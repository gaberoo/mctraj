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
    curState += trans.trans;
    curState.branches += trans.branchTrans;
    last_event_time += trans.time;
  }

  // =========================================================================

  int Trajectory::step(double maxTime, const void* pars, rng::RngStream* rng, 
                       bool noTree, bool adjZero) 
  {
    // get array of next event rates
    double nextTime = 0.0;
    double r = 0.0;

    int nextEvent = 0;
    int num_events = 0;

    double dw = 0.0;
    double rate = 0.0;
    double unrate = 0.0;

    int ret = 0;

    const TransitionType* nextTrans;

    // calculate all rates
    double totalRate = model->calculateTransRates(curState,transRates);

    while (ret <= 0) {
      if (totalRate <= 0.0) {
        time = maxTime;
        ret = 4;
        break;
      }

      // sample next event time
      rng->uniform(1,&r);
      nextTime =  -1.0/totalRate * log(r);

      if (time+nextTime < maxTime) {
        // sample next event type and calculate rate of the event
        nextEvent = model->chooseTransition(rng,transRates);
        rate = transRates[nextEvent] - ((nextEvent>0) ? transRates[nextEvent-1] : 0.0);

        // get the appropriate transition
        nextTrans = model->getTType(nextEvent);
        StateTransition st(nextTime,*nextTrans,curState,pars,nextEvent,time+nextTime);

        // get potential new state
        EpiState newState = curState;
        newState += st.trans;

        // check if new state is allowed
        dw = nextTrans->applyProb(newState,model->getPars());

        if (dw > 0.0 || ! adjZero) {
          // if yes,
          // add branch transformations
          if (! noTree) {
            curState.time = time + nextTime;
            nextTrans->applyBranch(curState,rng,st,pars);
            dw *= st.relprob;
          }

          // add transition to list
          if (store_trans) transitions.push_back(st);

          // adapt state
          addTransition(st);

          // multiply probability
          updateProb(dw);

          // advance time
          time += nextTime;
          ++num_events;

#ifdef DEBUG
          if (ret < 0) {
            cerr << "Legal event (" << nextEvent << ")."
                 << " w = " << dw << ", penalty = " << rate 
                 << ", total rate = " << totalRate << "." << endl;
          }
#endif

          ret = 1; // exit loop
        } else {
          // if no,
          // make this transition illegal in this particular step
#ifdef DEBUG
          cerr << "Illegal event (" << nextEvent << ")."
               << " w = " << dw << ", penalty = " << rate 
               << ", total rate = " << totalRate << "." << endl;
#endif
          // adjust cumulative rates
          for (size_t i = nextEvent; i < transRates.size(); ++i) {
            transRates[i] -= rate;
          }
          totalRate = transRates.back();

          // condition on not seeing this event
          unrate += rate;

          // redo loop with modified rates
          ret = -1;
        }
      } else {
#ifdef DEBUG
        if (ret < 0) {
          cerr << "No event. Total rate = " << totalRate << "." << endl;
        }
#endif

        nextTime = maxTime - time;
        time = maxTime;
        ret = 2; // exit loop
      }
    }
    
    if (totalRate < 0.0) {
#ifdef DEBUG
      cerr << "\033[1;31m";
      cerr << "Negative rate (" << ret << "): " << totalRate
           << "  ES = " << curState;
      cerr << "\033[0m" << endl;
#endif
      return -1;
    }

    // condition on the illegal events ('prob' is in logarithmic scale)
    if (unrate > 0.0) prob -= unrate*nextTime;

    return ret;
    // return num_events;
    // return (ret != 2) ? num_events : -2;
  }

  // =========================================================================

  double Trajectory::force(double nextTime, int nextEvent, 
                           const vector<int>& ids, rng::RngStream* rng,
                           const void* pars) 
  {
    // store this info so that it can be accessed by the model
    curState.curBranch = ids;

    // get the transition of the next event
    const TransitionType* tt = model->getObsType(model->mapType(nextEvent));

    // setup the transition
    StateTransition st(nextTime-lastEventTime(),*tt,getState(),pars,nextEvent);

    // probability that the event happened
    double dw = tt->applyRate(curState,pars);

    // apply the branch changes
    tt->applyBranch(curState,rng,st,pars);

    if (store_trans) transitions.push_back(st);

    // apply the state change
    addTransition(st);

    // update the time
    time = nextTime;

    return dw;
  }

  // =========================================================================

  int Trajectory::simulateTrajectory(double endTime, const void* pars, 
                                     rng::RngStream* rng, bool adjZero) 
  {
    int ret = 1;
    bool noTree = false;
    while (time < endTime && ret > 0) {
      // ret is the number of events in the step
      ret = step(endTime,pars,rng,noTree,adjZero);
    }
    return ret;
  }

  // =========================================================================

  double Trajectory::trajProb(const void* pars, int last) const {
    double w2;
    double dw = 1.0;

    EpiState x(curState);
    const TransitionType* tt;
    int ntrans = transitionCount();
    if (ntrans <= last) return 0.0;

    // Events that were not included in the tree
    if (ntrans > 0) {
      for (int i = ntrans-1; i > last; --i) {
        tt = getTrans(i).getType();
        w2 = tt->applyProb(x,pars);
        dw *= w2;
        x -= getTrans(i).getTrans();
//        if (w2 <= 0.0) {
//          cerr << getTrans(i).realTime() 
//               << "[" << i << "] >> " 
//               << tt->getName() << " " 
//               << tt->applyProb(curState,pars) << endl;
//        }
      }
    }

    return dw;
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

  void Trajectory::toTree(rng::RngStream* rng, vector<TreeNode>& tree,
                          const int lineageStates[]) const 
  {
    EpiState x(initialState);
    size_t ntrans = transitionCount();
    double time = 0.0;

    vector< vector<int> > states(model->n());

    size_t i, j;
    int branch_id = 0;
    int epi_id = 0;

    // Initialize
    for (i = 0; i < x.numStates(); ++i) {
      if (lineageStates[i]) {
        for (j = 0; j < (size_t) x[i]; ++j) {
          tree.push_back(TreeNode());
          branch_id = tree.size()-1;
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
          rng->uniform_int(1,&rand,0,states[2].size());

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
          rng->uniform_int(1,&rand,0,states[2].size());
          parent = states[2][rand];
          states[2].erase(states[2].begin()+rand);
          tree[parent].age = time;
          tree[parent].n_off = 0;

          p = st.getType()->applyProb(x,model->getPars());
          rng->uniform(1,&r);

          if (r < p) {
            add_extant(tree,parent);
            // tree[parent].extant_off = 2;
          }

          break;

        case 2:
          rng->uniform_int(1,&rand,0,states[1].size());
          parent = states[1][rand];
          states[1].erase(states[1].begin()+rand);
          tree[parent].age = time;
          tree[parent].n_off = 1;

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


