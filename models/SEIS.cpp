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
  es.branches.colWeight(w,1);

  // pick a branch for the infection to happen on
  int id = rng->pick(w.data(),w.size());

  if (id >= 0) {
#ifdef DEBUG
    cout << ascii::magenta << "    infBranch  "  << ascii::end
         << " on " << id << " :: " << es << " : "
         << es.branches.states[id].to_json() << endl;
#endif

    // add dummy transition
    st.addBranchTrans(id,1,1);
    st.lastBT().change.assign(2,0);
    st.lastBT().change[0] = 1;

  } else {
#ifdef DEBUG
    cout << ascii::magenta << "    infBranch  "  << ascii::end << " :: ";
    cout << "no branch was chosen!" << endl;
#endif
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
  int nI = es.branchState(es.cur(0),1);

#ifdef DEBUG
  cout << ascii::blue << "    infTreeObs" << ascii::end 
       << " :: " << nI << " total I in group." << endl;
#endif

  // get infectious force for the group
  trueRate = ep->beta*es[0]*nI/ep->N;
  // trueRate = (nI > 0) ? ep->beta*es[0]/ep->N : 0.0;

  return trueRate;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infBranchObs(const EpiState& es, rng::RngStream* rng, 
                                    StateTransition& st, const void* pars)
{
  // get number of infecteds in the branch group
  int nI = es.branchState(es.cur(0),1);

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

    st.addBranchTrans(es.cur(0),1,-1);
    st.lastBT().change.assign(2,0);
    st.lastBT().change[1] = -1;

    st.addBranchTrans(b1,es.branchCol(b1),1);
    st.lastBT().change.assign(2,0);
    st.lastBT().change[1] = 1;

    st.addBranchTrans(b2,es.branchCol(b2),0);
    st.lastBT().change.assign(2,0);
    st.lastBT().change[0] = 1;
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

/****************************************************************************/
/* RECOVERY FUNCTIONS *******************************************************/
/****************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars, double& trueRate) 
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

double MCTraj::SEISModel::recovTreeObs(const EpiState& es, const void* pars, double& trueRate) 
  // probability that a recovery was observed, i.e. a sampling
{
  EpiPars* ep = (EpiPars*) pars;

  // get number of I in group
  int nI = es.branchState(es.cur(0),1);

  // sampling rate for this group
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
  es.branches.colProb(w,1);

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
    int col = es.branches.getCol(id);

    if (es.branches.awake(id) && es.branches.all(id) == 1) {
      cerr << ascii::red;
      cerr << "Can't remove last remaining state from an awake branch!";
      cerr << ascii::end << endl;
      // st.relprob = -INFINITY;
    } else {
#ifdef DEBUG
    cout << ascii::magenta << "    recovBranch"  << ascii::end
         << " on " << id << " :: " << es << " : "
         << es.branches.states[id].to_json() << endl;
#endif

      st.addBranchTrans(id,col,col);
      st.lastBT().change.assign(2,0);
      st.lastBT().change[1] = -1;
    }

    // adjust for not picking the last branch
    st.relprob = 1.0 - adj/tot;
  } 
  else 
    // otherwise this means there was no possible recovery
  {
#ifdef DEBUG
    cout << ascii::red << "    recovBranch"  << ascii::end
         << " on " << id << " :: " << es << " : "
         << "couldn't choose a branch to recover!" << endl;
    cout << es.branches.to_json() << endl;
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
    if (es.branchState(es.cur(0),1) > 0) {
      st.addBranchTrans(es.cur(0),1,-1);
      st.lastBT().change.assign(2,0);
      st.lastBT().change[1] = -1;
    } else {
#ifdef DEBUG
      cout << ascii::red << "    recovBranchObs on " << es.cur(0) << " :: "
           << " not enough I in group: " << es.branchState(es.cur(0),1) 
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
  es.branches.colWeight(w,0);

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
       << es.branches.states[id].to_json() << endl;
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
  st.addBranchTrans(id,es.branchCol(id),1);
  st.lastBT().change.assign(2,0);
  st.lastBT().change[0] = -1;
  st.lastBT().change[1] = 1;

  if (es.branchCol(id) == 0) {
    st[3] = -1;
    st[4] = 1;
  }

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
      if (es.branchState(es.cur(0),0) > 0) 
        // there are exposed individuals in the branch cluster
      // if (es.branchCol(es.cur(0)) == 0)
      {
        st.addBranchTrans(es.cur(0),0,-1);

        // update branch states
        st.lastBT().change.assign(2,0);
        st.lastBT().change[0] = -1;
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
        st.lastBT().change.assign(2,0);
        st.lastBT().change[1] = 1;
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

double MCTraj::SEIS::sample_rho(const EpiState& es, rng::RngStream* rng, void* pars) const {
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



