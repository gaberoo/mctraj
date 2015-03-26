#include "Trajectory.h"

namespace MCTraj {
  Trajectory& Trajectory::copyState(const Trajectory& T) {
    time = T.time;
    last_event_time = T.time;
    initialState = T.curState;
    curState = T.curState;
    model = T.model;
    transRates = T.transRates;
    // transRates.resize(model->n());
    prob = T.prob;
    return *this;
  }

  // =========================================================================

  Trajectory& Trajectory::operator=(const Trajectory& T) {
    time = T.time;
    last_event_time = T.last_event_time;
    initialState = T.initialState;
    curState = T.curState;
    transitions = T.transitions;
    model = T.model;
    transRates = T.transRates;
    prob = T.prob;
    return *this;
  }

  // =========================================================================

  Trajectory& Trajectory::operator+=(const Trajectory& T) {
    transitions.insert(transitions.end(),
                       T.transitions.begin(),
                       T.transitions.end());
    return *this;
  }

  // =========================================================================

  void Trajectory::addTransition(const StateTransition& trans) {
    // transitions.push_back(trans);
    curState += trans.trans;
    curState.branches += trans.branchTrans;
    last_event_time += trans.time;
  }

  // =========================================================================

  int Trajectory::step(double maxTime, const void* pars, gsl_rng* rng, bool noTree) {
    // get array of next event rates
    double nextTime = 0.0;
    double totalRate = 0.0;
    double r = 0.0;
    int nextEvent = 0;
    int num_events = 0;

    // cerr << curState << " -- " << totalRate << endl;
    totalRate = model->calculateTransRates(curState,transRates);

    if (totalRate < 0.0) {
      return -1;
    } else if (totalRate == 0.0) {
      return 0;
    }

    r = gsl_rng_uniform(rng);
    nextTime =  -1.0/totalRate * log(r);
    nextEvent = model->chooseTransition(rng,transRates);

    if (time+nextTime < maxTime) {
      const TransitionType* nextTrans = model->getTType(nextEvent);
      // cerr << "S> " << nextTrans->getName() << endl;
      transitions.push_back(StateTransition(nextTime,*nextTrans,
                            curState,pars,nextEvent,time+nextTime));
      StateTransition& st = transitions.back();
//      StateTransition st(nextTime,*nextTrans,
//                         curState,pars,nextEvent,time+nextTime);
      if (! noTree) nextTrans->applyBranch(curState,rng,st,pars);
      addTransition(st);
      time += nextTime;
      ++num_events;
    } else {
      time = maxTime;
    }

    return num_events;
  }

  // =========================================================================

  double Trajectory::force(double nextTime, int nextEvent, 
                           const vector<int>& ids, 
                           gsl_rng* rng,
                           const void* pars) 
  {
    curState.curBranch = ids;

    const TransitionType* tt = model->getObsType(model->mapType(nextEvent));

    transitions.push_back(StateTransition(nextTime-lastEventTime(),*tt,getState(),pars,nextEvent));
    StateTransition& st = transitions.back();

    /* probability that the event happened */
    double dw = tt->applyRate(curState,pars);
    // cout << " ]]]]]]] " << nextEvent << " " << model->mapType(nextEvent) << " " << dw << endl;
    tt->applyBranch(curState,rng,st,pars);
    addTransition(st);
    time = nextTime;

    return dw;
  }

  // =========================================================================

  int Trajectory::simulateTrajectory(double endTime, const void* pars, 
                                     gsl_rng* rng) 
  {
    int ret = 1;
    // cerr << " SIM :: time " << " " << endTime << endl;
    while (time < endTime && ret > 0) {
      ret = step(endTime,pars,rng);
    }
    return ret;
  }

  // =========================================================================

  ostream& operator<<(ostream& out, const Trajectory& traj) {
    EpiState x(traj.initialState);
    double time(0.0);
    double dt(0.0);
//    int eventType(-1);
    out << setw(12) << time << " " << x << " " << endl;
    for (size_t i(0); i < traj.transitions.size(); ++i) {
      x += traj.transitions[i].getTrans();
      dt = traj.transitions[i].atTime();
      time += dt;
      out << setw(12) << time << " " << x << " ";
      out << " " << traj.transitions[i];
      out << endl;
    }
    return out;
  }

  // =========================================================================

  Trajectory::Trajectory(size_t n, istream* fh) {
    string line;
    getline(*fh,line);
    istringstream iss(line);
    iss >> time;
    EpiState cur(n);
    for (size_t i(0); i < n; ++i) iss >> cur[i];
    initialState = cur;
    EpiState last(cur);
    while (getline(*fh,line)) {
      transitions.push_back(StateTransition(n));
      for (size_t i(0); i < n; ++i) {
        iss >> cur[i];
        transitions.back()[i] = last[i]-cur[i];
      }
    }
    curState = cur;
  }

  // =========================================================================

  int Trajectory::lastType(int type) const {
    if (nTrans() > 0) {
      for (int i(nTrans()-1); i >= 0; --i) {
        if (transitions[i].etype() == type) {
          return i;
        }
      }
    }
    return -1;
  }

  // =========================================================================

  void Trajectory::printFromLast(size_t last) const {
    EpiState curState(getState());
    size_t ntrans = transitionCount();

    // Events that were not included in the tree
    if (ntrans > 0) {
      for (size_t i = ntrans-1; i > last; --i) {
        cout << curState << endl;
        StateTransition st = getTrans(i);
        curState -= st.getTrans();
      }
    }
  }

  // =========================================================================

  ostream& Trajectory::printFromFirst(ostream& out) const {
    EpiState x(initialState);
    size_t ntrans = transitionCount();
    double time = 0.0;

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

    writer.StartArray();

    writer.StartObject();
    writer.String("time");
    writer.Double(time);
    writer.String("state");
    writer.StartArray();
    for (size_t j = 0; j < x.numStates(); ++j) writer.Int(x[j]);
    writer.EndArray();
    writer.EndObject();

    for (size_t i = 0; i < ntrans; ++i) {
      x += transitions[i].getTrans();
      time += transitions[i].atTime();

      writer.StartObject();
      writer.String("time");
      writer.Double(time);
      writer.String("state");
      writer.StartArray();
      for (size_t j = 0; j < x.numStates(); ++j) writer.Int(x[j]);
      writer.EndArray();
      writer.EndObject();
    }

    writer.EndArray();

    out << buffer.GetString();

    return out;
  }

  // =========================================================================

  ostream& Trajectory::printBranches(ostream& out) const {
    EpiState x(initialState);

    double time(0.0);
    double dt(0.0);

    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);

    writer.StartArray();
    for (size_t i(0); i < transitions.size(); ++i) {
      x += transitions[i].getTrans();
      dt = transitions[i].atTime();
      time += dt;
      for (size_t j(0); j < transitions[i].branchTrans.size(); ++j) {
        writer.StartObject();

        writer.String("time");
        writer.Double(time);

        writer.String("id");
        writer.Int(transitions[i].branchTrans.at(j).id);

        writer.String("old");
        writer.Int(transitions[i].branchTrans.at(j).old_color);

        writer.String("new");
        writer.Int(transitions[i].branchTrans.at(j).new_color);

        writer.EndObject();
      }
    }
    writer.EndArray();
    out << buffer.GetString();
    return out;
  }

  // =========================================================================

  ostream& Trajectory::printTxt(ostream& out) const {
    EpiState x(initialState);
    size_t ntrans = transitionCount();
    double time = 0.0;

    out << setw(12) << time << x << endl;
    for (size_t i = 0; i < ntrans; ++i) {
      x += transitions[i].getTrans();
      time += transitions[i].atTime();
      out << setw(12) << time << " "
          << x << " "
          << transitions[i] << " "
          <<endl;
    }

    return out;
  }

  // =========================================================================

  void Trajectory::toTree(gsl_rng* rng, vector<TreeNode>& tree) const {
    EpiState x(initialState);
    size_t ntrans = transitionCount();
    double time = 0.0;

    vector< vector<int> > states(model->n());

    size_t i, j;
    int branch_id = 0;
    int epi_id = 0;

    // Initialize
    int lineageStates[] = { 0, 1, 1, 0, 0 };

    for (i = 0; i < x.numStates(); ++i) {
      if (lineageStates[i]) {
        for (j = 0; j < (size_t) x[i]; ++j) {
          tree.push_back(TreeNode());
          branch_id = tree.size()-1;
          tree.back().code = branch_id;
          tree.back().code = branch_id;
          tree.back().epi_id = epi_id++;
          tree.back().epi_state = i;
          states[i].push_back(branch_id);
        }
      }
    }

    int rand;
    int parent;

    double p = 1.0;
    double r;
    StateTransition st;

    for (size_t i = 0; i < ntrans; ++i) {
      st = transitions[i]; 
      x += st.getTrans();
      time += st.atTime();

      switch (st.etype()) {
        case 1:
          // choose a parent
          rand = gsl_rng_uniform_int(rng,states[2].size());

          parent = states[2][rand];
          states[2].erase(states[2].begin()+rand);

          tree[parent].age = time;
          tree[parent].n_off = 2;

          // old infected
          tree.push_back(TreeNode());
          branch_id = tree.size()-1;

          tree[branch_id].code = branch_id;
          tree[branch_id].parent = tree[parent].code;
          tree[branch_id].age = -1.0;
          tree[branch_id].ancestor_age = tree[parent].age;
          tree[branch_id].epi_id = tree[parent].epi_id;
          tree[branch_id].epi_state = tree[parent].epi_state;

          states[2].push_back(branch_id);
          tree[parent].off.push_back(branch_id);

          // new infected
          tree.push_back(TreeNode());
          branch_id = tree.size()-1;

          tree[branch_id].code = branch_id;
          tree[branch_id].parent = tree[parent].code;
          tree[branch_id].age = -1.0;
          tree[branch_id].ancestor_age = tree[parent].age;
          tree[branch_id].epi_id = epi_id++;
          tree[branch_id].epi_state = 1;

          states[1].push_back(branch_id);
          tree[parent].off.push_back(branch_id);

          break;

        case 0:
          rand = gsl_rng_uniform_int(rng,states[2].size());
          parent = states[2][rand];
          states[2].erase(states[2].begin()+rand);
          tree[parent].age = time;

          p = st.getType()->applyProb(x,model->getPars());
          r = gsl_rng_uniform(rng);

          if (r < p) tree[parent].extant_off = 2;

          break;

        case 2:
          rand = gsl_rng_uniform_int(rng,states[1].size());
          parent = states[1][rand];
          states[1].erase(states[1].begin()+rand);
          tree[parent].age = time;

          tree.push_back(TreeNode());
          branch_id = tree.size()-1;
          states[2].push_back(branch_id);
          tree[branch_id].code = branch_id;
          tree[branch_id].parent = tree[parent].code;
          tree[branch_id].age = -1.0;
          tree[branch_id].ancestor_age = tree[parent].age;
          tree[branch_id].epi_id = tree[parent].epi_id;
          tree[branch_id].epi_state = 2;
          tree[parent].off.push_back(branch_id);

          break;
      }
    }

    // print the remaining branches
    for (i = 0; i < x.numStates(); ++i) {
      if (lineageStates[i]) {
        for (j = 0; j < states[i].size(); ++j) {
          tree[states[i][j]].age = time;
        }
      }
    }
  }
}


