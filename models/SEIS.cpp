#include "SEIS.h"

double MCTraj::SEISModel::infRateFun(const EpiState& es, const void* pars, double& trueRate) 
  // The rate at which any new infections happen.
{
  EpiPars ep = *(EpiPars*) pars;
  trueRate = (es[2] > 0) ? ep.beta*es[0]*es[2]/ep.N : 0.0;
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::infTreeProb(const EpiState& es, const void* pars) 
  // The probability that a new infection was not visible in the tree.
{
  return 1.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infValidate(const EpiState& es, const void* pars) {
  return 1;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infBranch(const EpiState& es, rng::RngStream* rng, 
                                 StateTransition& st, const void* pars)
{
  // get the total number of infecteds from the branches 
  // (this should be equal to I)
  vector<double> w;
  es.branches.colWeight(w,2);

  // pick a branch for the infection to happen on
  int id = rng->pick(w.data(),w.size());

  if (id >= 0) {
#ifdef DEBUG
    cout << ascii::magenta << "    infBranch  "  << ascii::end
         << " on " << id << " :: " << es << " : "
         << es.branches.state_to_json(id) << endl;
#endif

    // add dummy transition
    st.addBranchTrans(id,es.branchCol(id),es.branchCol(id));
    st.lastBT().change.assign(3,0);
    st.lastBT().change[0] = 1;  // overall counter
    st.lastBT().change[1] = 1;  // color counter

  } else {
#ifdef DEBUG
    cout << ascii::magenta << "    infBranch  "  << ascii::end << " :: ";
    cout << "no branch was chosen!" << endl;
#endif
    st.relprob = -INFINITY;
  }

  return 0;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::infTreeObs(const EpiState& es, const void* pars, double& trueRate) 
  // The probability that an observed infection event happens given the
  // current population state
{
  EpiPars* ep = (EpiPars*) pars;

  // get number of infecteds in the branch group
  double nE = es.branchState(es.cur(0),1);
  double nI = es.branchState(es.cur(0),2);

  // get infectious force for the group:
  //    rate of infection on a single branch
  //  x probability that the branch is in state I
  // trueRate = ep->beta*es[0]/ep->N * nI/(nI+nE);
  trueRate = ep->beta*es[0]/ep->N * nI;

#ifdef DEBUG
  cout << ascii::blue << "    infTreeObs" << ascii::end 
       << " :: " << nI << " total I in group. "
       << "Rate = " << setw(8) << trueRate << " | " << es << " | "
       << es.branches.state_to_json(es.cur(0))
       << endl;
#endif

  return trueRate;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infBranchObs(const EpiState& es, rng::RngStream* rng, 
                                    StateTransition& st, const void* pars)
{
  // get number of infecteds in the branch group
  int nI = es.branchState(es.cur(0),2);

  int b1 = es.cur(1);
  int b2 = es.cur(2);

  double r;
  rng->uniform(1,&r);

  if (r < 0.5) { 
    int tmp = b1; b1 = b2; b2 = tmp;
  }

  if (nI > 0) {
#ifdef DEBUG
    cout << "        infBranchObs ("
         << es.cur(0) << "," << b1 << "," << b2 << ")" << endl;
#endif

    st.addBranchTrans(es.cur(0),es.branchCol(es.cur(0)),-1);
    st.lastBT().change.assign(3,0);
    st.lastBT().change[2] = -1;

    st.addBranchTrans(b1,es.branchCol(b1),1);
    st.lastBT().change.assign(3,0);
    st.lastBT().change[0] = 1;
    st.lastBT().change[2] = 1;

    st.addBranchTrans(b2,es.branchCol(b2),0);
    st.lastBT().change.assign(3,0);
    st.lastBT().change[0] = 1;
    st.lastBT().change[1] = 1;
  } else {
#ifdef DEBUG
    cout << ascii::red << "        infBranchObs" << ascii::end 
         << " :: " << "no I in group!" << endl;
#endif

    st.relprob = -INFINITY;
    return -1;
  }

  return 0;
}

/****************************************************************************\
|* RECOVERY FUNCTIONS *******************************************************|
\****************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars, 
                                       double& trueRate) 
  // rate at which recovery occurs
{
  EpiPars* ep = (EpiPars*) pars;
  trueRate = (ep->mu+ep->psi)*es[2];
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeProb(const EpiState& es, const void* pars) 
  // the probability that the event did not result in a sampling
{
  EpiPars* ep = (EpiPars*) pars;
  return ep->mu/(ep->mu+ep->psi);
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::recovValidate(const EpiState& es, const void* pars) 
{
  return 1;
}

//---------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeObs(const EpiState& es, const void* pars, 
                                       double& trueRate) 
  // probability that a recovery was observed, i.e. a sampling
{
  EpiPars* ep = (EpiPars*) pars;

  // get number of I in group
  double nE = es.branchState(es.cur(0),1);
  double nI = es.branchState(es.cur(0),2);

  // sampling rate for this group
  // trueRate = ep->psi * nI/(nI+nE);
  trueRate = ep->psi * nI;

#ifdef DEBUG
  cout << ascii::blue << "    recovTreeObs" << ascii::end 
       << " :: " << nI << " total I in group." << endl;
#endif

  return trueRate;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::recovBranch(const EpiState& es, rng::RngStream* rng, 
                                   StateTransition& st, const void* pars)
{
  // get number of I in all branches
  vector<double> w;
  es.branches.colProb(w,2);

  // check for last branch in live group
  double tot = 0.0;
  double adj = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    tot += w[i];
    if (es.branchAwake(i) && es.branchSize(i) == 1) {
      // this is the last branch in the group
      adj += w[i];
      w[i] = 0.0;
    }
  }
  // make cumulative sum
  for (size_t i = 1; i < w.size(); ++i) w[i] += w[i-1];

  int id = rng->pick(w.data(),w.size());
  if (id >= 0) 
    // make sure a valid branch was found
  {
    if (es.branchAwake(id) && es.branches.all(id) == 1) {
      cerr << ascii::red;
      cerr << "Can't remove last remaining state from an awake branch!";
      cerr << ascii::end << endl;
      // st.relprob = -INFINITY;
    } else {
#ifdef DEBUG
    cout << ascii::magenta << "    recovBranch"  << ascii::end
         << " on " << id << " :: " << es << " : "
         << es.branches.state_to_json(id) << endl;
#endif

      st.addBranchTrans(id,es.branchCol(id),es.branchCol(id));
      st.lastBT().change.assign(3,0);
      st.lastBT().change[2] = -1;
    }

    // adjust for not picking the last branch
    st.relprob = log(1.0-adj/tot);
  } 
  else 
    // otherwise this means there was no possible recovery
  {
#ifdef DEBUG
    cout << ascii::red << "    recovBranch"  << ascii::end
         << " on " << id << " :: " << es << " : "
         << "couldn't choose a branch to recover!" << endl;
    cout << "      weight array" << endl;
    cout << "        ";
    for (size_t i = 0; i < w.size(); ++i) {
      cout << setw(3) << i << " ";
    }
    cout << endl;
    cout << "        ";
    for (size_t i = 0; i < w.size(); ++i) {
      cout << setw(3) << rng::dprob(i,w.data()) << " ";
    }
    cout << endl;
#endif

    st.relprob = -INFINITY;
    return -1;
  }

  return 0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::recovBranchObs(const EpiState& es, rng::RngStream* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
    return -1;
  } else {
    if (es.branchState(es.cur(0),2) > 0) {
      st.addBranchTrans(es.cur(0),es.branchCol(es.cur(0)),-1);
      st.lastBT().change.assign(3,0);
      st.lastBT().change[2] = -1;
    } else {
#ifdef DEBUG
      cout << ascii::red << "    recovBranchObs on " << es.cur(0) << " :: "
           << " not enough I in group: " << es.branchState(es.cur(0),2) 
           << ascii::end << endl;
#endif
      return -2;
    }
  }
  return 0;
}

/************************************************************************/
/* TRANSITIONS **********************************************************/
/************************************************************************/

double MCTraj::SEISModel::transRateFun(const EpiState& es, const void* pars, double& trueRate) 
  // transition rate in the simulation
{
  EpiPars* ep = (EpiPars*) pars;
  trueRate = ep->gamma * es[1];
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::transTreeProb(const EpiState& es, const void* pars) 
  // probability that the transition didn't happen on the tree
{
  return 1.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::transValidate(const EpiState& es, const void* pars) 
  // some validation function I can't remember
{
  return 1;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::transTreeObs(const EpiState& es, const void* pars, double& trueRate) 
  // probability that the transition was observed, given that it's on the tree
{
  trueRate = 1.0;
  return 1.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::transBranch(const EpiState& es, rng::RngStream* rng, 
                                   StateTransition& st, const void* pars)
  // decide whether the transition happened on a branch or not, and if yes,
  // apply that transition to the branch
{
  // get the weight of each branch
  vector<double> w;
  es.branches.colWeight(w,1);

  // pick branch for the transition to happen on
  int id = rng->pick(w.data(),w.size());

  if (id < 0) {
#ifdef DEBUG
    cout << "Couldn't choose a branch to transition on!!! kE = 0?" << endl;
#endif
    return 0;
  }

#ifdef DEBUG
  cout << ascii::magenta << "    transBranch" << ascii::end
       << " on " << id << " :: " << es << " : "
       << es.branches.state_to_json(id) << endl;
//  cout << "      weight array" << endl;
//  cout << "        ";
//  for (size_t i = 0; i < w.size(); ++i) {
//    cout << setw(3) << i << " ";
//  }
//  cout << endl;
//  cout << "        ";
//  for (size_t i = 0; i < w.size(); ++i) {
//    cout << setw(3) << rng::dprob(i,w.data()) << " ";
//  }
//  cout << endl;
#endif

  // add dummy transition
  st.addBranchTrans(id,es.branchCol(id),es.branchCol(id));
  st.lastBT().change.assign(3,0);
  st.lastBT().change[1] = -1;
  st.lastBT().change[2] = 1;

  return 0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::transBranchObs(const EpiState& es, rng::RngStream* rng, 
                                      StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
    // check for node A information
    if (es.cur(0) >= 0) {
      // check the state of the current branch
      if (es.branchState(es.cur(0),1) > 0) 
        // there are exposed individuals in the branch cluster
      // if (es.branchCol(es.cur(0)) == 0)
      {
        st.addBranchTrans(es.cur(0),0,-1);

        // update branch states
        st.lastBT().change.assign(3,0);
        st.lastBT().change[1] = -1;
#ifdef DEBUG
        // cout << "transBranchObs :: Putting branch " << es.cur(0) << " to sleep." << endl;
#endif
      } else {
#ifdef DEBUG
        cout << "transBranchObs :: Node A doesn't have the correct color (" 
             << es.cur(0) << ","
             << es.cur(1) << ","
             << es.cur(2) << ") : " 
             << es.branchCol(es.cur(0)) << endl;
#endif
        return -1;
      }
    }

    // check for node B information
    if (es.cur(1) >= 0) {
      if (es.branchCol(es.cur(1)) == -1) 
        // new branch is asleep
      {
        st.addBranchTrans(es.cur(1),-1,1);

        // update branch states
        st.lastBT().change.assign(3,0);
        st.lastBT().change[2] = 1;
#ifdef DEBUG
        // cout << "transBranchObs :: Waking branch " << es.cur(1) << "." << endl;
#endif
      } else {
#ifdef DEBUG
        cout << "transBranchObs :: Node B isn't asleep (" 
             << es.cur(0) << ","
             << es.cur(1) << ","
             << es.cur(2) << ") : " 
             << es.branchCol(es.cur(0)) << endl;
#endif
        return -2;
      }
    }
    // cerr << st.to_json() << endl;
  }
  return 0;
}

/************************************************************************/

double MCTraj::SEIS::sample_rho(const EpiState& es, rng::RngStream* rng, 
                                void* pars) const 
{
  int k = es[3]+es[4];
  int I = es[1]+es[2];
  double w = 0.0;
  if (rho >= 1.0 && I == k) {
    w = 1.0;
  } else if (rho < 1.0 && I >= k) {
    w = gsl_ran_binomial_pdf(k,rho,I);
  }
  return w;
}

/************************************************************************/

bool MCTraj::SEIS::validState(const EpiState& es) const {
  if (es[1] < es[2]) return false;
  if (es[2] < es[4]) return false;
  if (es[0] < 0 || es[1] < 0 || es[2] < 0 || es[3] < 0 || es[4] < 0) return false;
  return true;
}

/************************************************************************/

double MCTraj::SEISModel::calcBranchPotentials
  (const EpiState& es, const void* pars, int color, double* rates)
// potentials on the branches
{
  EpiPars* ep = (EpiPars*) pars;

  size_t n = es.branches.nAlive();
  if (n == 0) return 0.0;

  double totRate = 0.0;
  double rate = 0.0;
  size_t i = 0;
  double btime = 0.0;
  double dt = 0.0;

  int id = es.branches.alive[0];

  // if the branch is in state E, then transition is possible with
  // rate as function of the distance from the next internal node
  if (es.branches.getCol(id) == color) {
    dt = ep->tree->branches.at(id).time - es.time;
    rate = exp(1.0/(dt + ep->alpha));
  } else {
//#ifdef DEBUG
//    cout << "wrong color: " << i << " " << id << " " 
//         << es.branches.getCol(id) << " " << color << endl;
//#endif
    rate = 0.0;
  }

  if (rates != NULL) rates[i] = rate;

  for (i = 1; i < n; ++i) {
    id = es.branches.alive[i];
    if (es.branches.getCol(id) == color) {
      btime = ep->tree->branches.at(id).time;
      dt = ep->tree->branches.at(id).time - es.time;
      if (btime <= es.time) {
        // this should never occur !
        cerr << "next branch time on " << id << " is before current time!" << endl
             << "   --> " << btime << " <> " << es.time << endl;
      }
      rate = exp(1.0/(dt + ep->alpha));
      if (rate == INFINITY) {
        cerr << "Rate is infinite!" << endl
             << "   " << i << " " << id << " " << btime << " " << es.time 
             << " " << exp(1.0/(dt + ep->alpha)) << endl;
      }
      if (rates != NULL) rates[i] = rate + rates[i-1];
      totRate += rate;
    } else {
//#ifdef DEBUG
//      cout << "wrong color: " << i << " " << id << " " 
//           << es.branches.getCol(id) << " " << color << endl;
//#endif
      if (rates != NULL) rates[i] = rates[i-1];
    }
  }

  return totRate;
}

/************************************************************************/

void MCTraj::SEIS::toTree(const Trajectory& traj,
                          rng::RngStream* rng, 
                          vector<TreeNode>& tree) const 
{
  EpiState x(traj.initState());
  size_t ntrans = traj.transitionCount();
  double time = 0.0;

  const Model* model = traj.getModel();
  vector< vector<int> > states(model->n());

  size_t i, j;
  int branch_id = 0;
  int epi_id = 0;

  int lineageStates[] = { 0, 1, 1, 0, 0 };

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
    st = traj.getTrans(i); 
    x += st.getTrans();
    time += st.atTime();

    cout << st.to_json() << endl;
    cout << "Next state = " << x.to_json() << endl;

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

