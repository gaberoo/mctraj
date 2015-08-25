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
  double E = es[1];
  double I = es[2];
  double kE = es[3];
  double kI = es[4];

  double p = 1.0 - kI/I * kE/E;
  // double p = 1.0;

#ifdef DEBUG
  cerr << "    treeProbInf   :: ES = " << es << " | p = " << p << endl;
#endif
  return (I >= kI && E >= kE) ? p : 0.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infValidate(const EpiState& es, const void* pars) {
  return 1;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::infTreeObs(const EpiState& es, const void* pars, double& trueRate) 
  // The probability that an observed infection event happens given the
  // current population state
{
  EpiPars ep = *(EpiPars*) pars;

  double S = es[0];
  // double E = es[1];
  double I = es[2];

  double lambda = ep.beta * S /ep.N;
  // double lambda = ep.beta;

  int branch_color = 0;

  double x = 0.0;

  // make sure there are any infecteds
  if (I <= 0.0) {
#ifdef DEBUG
    cerr << ascii::b_red 
         << "    treeObsInf    :: ES = " << es << " | no infecteds!" 
         << ascii::end << endl;
#endif
    return 0.0;
  }

  if (es.curBranch.size() > 0) 
    // check if branching information was specified
  {
    // get the color of the specified branch
    branch_color = es.branches.getCol(es.curBranch[0]);

    // if the branch is in state I, then infection is possible
    x = (branch_color == 1) ? lambda : 0.0;
    // x = (branch_color == 1) ? lambda/(E+1.) : 0.0;

#ifdef DEBUG
    if (x <= 0.0) {
      cerr << ascii::b_red
           << "    treeObsInf    :: ES = " << es 
           << " | " << es.curBranch[0] << " -- " << branch_color;
      cerr << ascii::end << endl;
      for (size_t j = 0; j < es.branches.alive.size(); ++j) {
        cerr << "        " << j << " " 
             << es.branches.alive[j] << " "
             << es.branches.colors[es.branches.alive[j]] << endl;
      }
      cerr << ascii::end << endl;
    }
#endif

    trueRate = x;
    return x;
  } else {

#ifdef DEBUG
    cerr << "    treeObsInf    :: no node information! ES = " << es << endl;
#endif

    trueRate = 1.0;
    return 1.0;
  }
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infBranch(const EpiState& es, rng::RngStream* rng, 
                                 StateTransition& st, const void* pars)
{
  EpiPars* ep = (EpiPars*) pars;

  // double E = es[1];
  double I = es[2];
  // double kE = es[3];
  double kI = es[4];

  double r; 
  rng->uniform(1,&r);

  // the probability of picking a tree branch is kI/I
  if (r < kI/I) {
    double p = 1.0;
    double dt = INFINITY;

    // get a random branch in state I (I=1)
    int id = es.branches.random_color(rng,1);

    if (id >= 0) {
      if (ep->tree != NULL) {
#ifdef DEBUG
        cerr << "        infBranch :: " << id << " (" << es.branches.getCol(id) << ")" << endl;
//        cerr << id << "/" << ep->tree->branches.size() << endl;
//        cerr << "Infection happened on branch "
//            << ep->tree->branches[id].id << ". " << endl;
//        cerr << "   Time left before next event = "
//            << ep->tree->branches[id].time - es.time << endl;
//        cerr << "   Type of next event = "
//            << ep->tree->branches[id].type << endl;
#endif

        // the probability that the new individual becomes this 
        // branch is 0.5 * (1 + exp(-dt)), where dt is the time 
        // to the next interval
        dt = ep->tree->branches[id].time - (es.time+es.nextTime);

        switch (ep->tree->branches[id].type) {
          case 0:
          case 1:
          case 2: p = 1.0 + exp(-dt); break;
          default: p = 1.0; break;
        }

#ifdef DEBUG
        cerr << "        Choice bias: p = " << 0.5*p << endl;
#endif

      }

      rng->uniform(1,&r);
      if (r < 0.5*p) {
        st.relprob = -log(p);
#ifdef DEBUG
        cerr << "        Branch remains in state (r = " << r << ")" << endl;
#endif
      } else {
        // the color change happens on the branch
        st.branchTrans.push_back(BranchStateChange(id,1,0));
        st[3] = 1;
        st[4] = -1;
        st.relprob = -log(2.0-p);
#ifdef DEBUG
        cerr << "        Branch changes state (r = " << r << ")" << endl;
#endif
      }
    } else {
      // otherwise raise an error
      cerr << "     [I] No branch of color 'I' exists! | ES = " << es << endl;
      cerr << es.branches.to_json() << endl;
      return -1;
    }
  }
  return 0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infBranchObs(const EpiState& es, rng::RngStream* rng, 
                                    StateTransition& st, const void* pars)
{
  EpiPars* ep = (EpiPars*) pars;

  // choose which branch becomes E and I

  int b0 = es.curBranch.at(0);
  int b1 = es.curBranch.at(1);
  int b2 = es.curBranch.at(2);

  if (es.branches.getCol(b0) != 1) 
    // the branch has the wrong color
  {
#ifdef DEBUG
    cerr << ascii::red 
         << "     Branch " << b0 << " is not in state I: " 
         << es.branches.getCol(b0) 
         << ascii::end << endl;
#endif
    return -1;
  }

  double p = 0.5;

  // Choose which branch keeps the infected. Priority is given when the next
  // event on a branch requires and infected class
  {
    double dt1 = ep->tree->branches[b1].time - (es.time+es.nextTime);
    double p1 = 1.0;
    switch (ep->tree->branches[b1].type) {
      case 0: // p1 = 1-exp(-dt1); break;
      case 1: p1 = 1+exp(-dt1); break;
      default: p1 = 1.0; break;
    }

    double dt2 = ep->tree->branches[b2].time - (es.time+es.nextTime);
    double p2 = 1.0;
    switch (ep->tree->branches[b2].type) {
      case 0: // p2 = 1-exp(-dt2); break;
      case 1: p2 = 1+exp(-dt2); break;
      default: p2 = 1.0; break;
    }

    p = p1/(p1+p2);
  }

  // int a, b;
  int tmp;
  double r;
  rng->uniform(1,&r);
  if (r < p) { 
    // if (ep->alpha >= 0.0) 
    { st.relprob = -M_LN2 - log(p); }
  } else { 
    tmp = b1; b1 = b2; b2 = tmp; // swap branches
    // if (ep->alpha >= 0.0) 
    { st.relprob = -M_LN2 - log(1.0-p); }
  }

  // get colors of the branches
  int col1 = es.branches.getCol(b1);
  int col2 = es.branches.getCol(b2);

#ifdef DEBUG
  cerr << "        infBranchObs :: (" << b0 << "," << b1 << "," << b2 << ")" << endl;
#endif

  if (es.curBranch.size() < 3) {
    cerr << "Need node information to color tree!" << endl;
  } else {
    if (es.curBranch[2] < 0) {
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(1),-1,1));
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
      st[1] = 0;  /* E ->   */
      st[2] = 1;  /*   -> I */
      st[3] = 0;
      st[4] = 1;  /* kI + 1 */
    } else {
      // cerr << es.curBranch.at(a) << " >> "
      //      << es.curBranch.at(b) << " || "
      //      << es.curBranch.at(0) << endl;
      st.branchTrans.push_back(BranchStateChange(b1,col1,1));
      st.branchTrans.push_back(BranchStateChange(b2,col2,0));
      st.branchTrans.push_back(BranchStateChange(b0,1,-1));
    }
  }

  return 0;
}

/****************************************************************************/
/* RECOVERY FUNCTIONS *******************************************************/
/****************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars, double& trueRate) 
{
  EpiPars ep = *(EpiPars*) pars;

  // only allow recovery in the non-tree lineages
  // return (ep.mu+ep.psi)*(es[2]-es[4]);

  trueRate = (ep.mu+ep.psi)*es[2];
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeProb(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.mu);

  double I  = es[2]+1.0;
  double kI = es[4];

  double dw = (I > kI) ? (1.-s)*(1.-kI/I) : 0.0;
  // double dw = (I > kI) ? (1.-s) : 0.0;

  if (es[2] < es[4] || es[1] < es[3]) dw = 0.0;

#ifdef DEBUG
  if (dw <= 0.0) cerr << ascii::red;
  cerr << "     treeProbRecov :: ES = " << es << " | w = " << dw 
       << " >> " << I << " --> " << kI 
       << ascii::end << endl;
#endif

  return dw;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::recovValidate(const EpiState& es, const void* pars) {
  return 1;
}

//---------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeObs(const EpiState& es, const void* pars, double& trueRate) 
{
  EpiPars ep = *(EpiPars*) pars;

  int id = es.curBranch[0];
  trueRate = 0.0;

  if (es.curBranch.size() > 0) 
    // check for branching information
  {
    if (es.branches.getCol(id) == 1) {
      trueRate = ep.psi*es[2];
      // trueRate = ep.psi*es[4];
      // trueRate = ep.psi; // <-- is this really just 'psi' ?
    } else {
#ifdef DEBUG
      cerr << "     treeObsRecov :: ES = " << es 
           << " | " << id << " @ " << es.branches.getCol(id) << endl;
#endif
    }
  } else {
    cerr << "     treeObsRecov :: no node information! ES = " << es << endl;
  }

  return trueRate;
}

/************************************************************************/

int MCTraj::SEISModel::recovBranchObs(const EpiState& es, rng::RngStream* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
    return -1;
  } else {
    if (es.branches.getCol(es.curBranch.at(0)) == 1) {
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
    } else {
#ifdef DEBUG
      cerr << "     Branch " << es.curBranch.at(0) << " is wrong color for (I->0): " 
           << es.branches.getCol(es.curBranch.at(0)) << endl;
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
  double E = es[1];

  trueRate = ep->gamma * E;

  if (ep->alpha >= 0.0) {
    double kE = es[3];
    double branchRates = 0.0;
    // cerr << "Starting CBP..." << flush;
    if (kE > 0) branchRates = calcBranchPotentials(es,pars,0);
    // cerr << "done. bR = " << branchRates << endl;
    return ep->gamma * ((E-kE) + branchRates);
  } else {
    return trueRate;
  }
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::transTreeProb(const EpiState& es, const void* pars) 
// probability that the transition didn't happen on the tree
{
  double E = es[1];
  double kE = es[3];
  // return (E >= kE) ? (E-kE)/E : 0.0;
  // return (E >= kE) ? 1.0 : 0.0;
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
  // debug("treeObsTrans :: ES = (%d,%d,%d,%d,%d)",es[0],es[1],es[2],es[3],es[4]);
  trueRate = 1.0;
  return 1.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::transBranch(const EpiState& es, rng::RngStream* rng, 
                                   StateTransition& st, const void* pars)
// decide whether the transition happened on a branch or not, and if yes,
// apply that transition to the branch
{
  EpiPars* ep = (EpiPars*) pars;
  // BranchStates braches(es.branches);

  double E = es[1];
  double kE = es[3];

  int n = 0;  // alive branches
  int i = 0;  // chosen rate

  if (ep->alpha >= 0.0) 
  {
    // cerr << "Using transition potentials..." << endl;

    double* rates = NULL;

    if (kE > 0) {
      n = es.branches.alive.size();
      rates = new double[n+1];
      calcBranchPotentials(es,pars,0,rates);
      rates[n] = E-kE + rates[n-1];
      i = rng->pick(rates,n+1);

#ifdef DEBUG
      if (i < 0) {
        cerr << "Can't pick a transition!" << endl;
        int id = 0;
        for (size_t j = 0; j < n; ++j) {
          id = es.branches.alive[j];
          cerr << "   " 
               << id << " < "
               << setw(7) << ep->tree->branches[id].time 
               << " => " << es.branches.colors[id] 
               << " >> " << j << " " << rates[j] << endl;
        }
        cerr << "   " 
             << 0 << " < "
             << setw(7) << 0.0 
             << " => " <<  -1
             << " >> " << n << " " << rates[n] << endl;
      }
      cerr << "[" << es.time << "/" << es.time + es.nextTime 
           << "] Transition on branch " << i << " of " << n
           << " (E=" << E << ",kE=" << kE << ")" 
           << endl;
#endif

      double logp = ep->gamma*es.nextTime*(rates[n-1]-kE);
      double rate = (i < n) ? (rates[i] - ((i>0) ? rates[i-1] : 0.0)) : 1.0;
      st.relprob = logp - log(rate);

#ifdef DEBUG
      if (logp-log(rate) > 0.0) {
        cerr << "simulated trajectory less likely!" << endl;
        cerr << i << "/" << n << " " 
             << setw(8) << es.time << " " 
             << setw(8) << es.nextTime << " " 
             << setw(8) << rate << " " 
             << setw(12) << logp << " "
             << endl;
        for (size_t j = 0; j < n; ++j) {
          int id = es.branches.alive[j];
          cerr << "   " 
               << id << " < "
               << setw(7) << ep->tree->branches[id].time 
               << " =>  " << es.branches.colors[id] 
               << " >> " << j << " " << rates[j] << endl;
        }
        cerr << "   " 
             << 0 << " < "
             << setw(7) << 0.0 
             << " => " <<  -1
             << " >> " << n << " " << rates[n] << endl;
      }
#endif

      if (i < n) 
        // check if the transition happened on a brach, i.e. no in the
        // last compartment in the rates
      {
        // int id = es.branches.random_color(rng,0);
        // cerr << "Trans on branch. " << r << "/" << id << " " << es << endl;
        int id = es.branches.alive.at(i);
        if (id >= 0) {
          st.branchTrans.push_back(BranchStateChange(id,0,1));
          st[3] = -1;
          st[4] = 1;
        } else {
          // cerr << "[T] Color doesn't exist (" << id << ") !" << endl;
          // cerr << es.branches.to_json() << endl;
          return -1;
        }
      }

      delete[] rates;
    }
  } else {
    double r = 0.0;
    rng->uniform(1,&r);
    if (r < kE/E) {
      int id = es.branches.random_color(rng,0);
      // cerr << "Trans on branch. " << r << "/" << id << " " << es << endl;
      if (id >= 0) {
        st.branchTrans.push_back(BranchStateChange(id,0,1));
        st[3] = -1;
        st[4] = 1;
      } else {
        // cerr << "[T] Color doesn't exist (" << id << ") !" << endl;
        // cerr << es.branches.to_json() << endl;
        return -1;
      }
    }
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
    if (es.curBranch.at(0) >= 0) {
      if (es.branches.getCol(es.curBranch.at(0)) == 0) {
        st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),0,-1));
      } else {
#ifdef DEBUG
        cerr << "transBranchObs :: Node A doesn't have the correct color (" 
             << es.curBranch.at(0) << ","
             << es.curBranch.at(1) << ","
             << es.curBranch.at(2) << ") : " 
             << es.branches.getCol(es.curBranch.at(0)) << endl;
#endif
        return -1;
      }
    }

    // check for node B information
    if (es.curBranch.at(1) >= 0) {
      if (es.branches.getCol(es.curBranch.at(1)) == -1) {
        st.branchTrans.push_back(BranchStateChange(es.curBranch.at(1),-1,1));
      } else {
#ifdef DEBUG
        cerr << "transBranchObs :: Node B isn't asleep (" 
             << es.curBranch.at(0) << ","
             << es.curBranch.at(1) << ","
             << es.curBranch.at(2) << ") : " 
             << es.branches.getCol(es.curBranch.at(0)) << endl;
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

  size_t n = es.branches.alive.size();
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
//    cerr << "wrong color: " << i << " " << id << " " 
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
//      cerr << "wrong color: " << i << " " << id << " " 
//           << es.branches.getCol(id) << " " << color << endl;
//#endif
      if (rates != NULL) rates[i] = rates[i-1];
    }
  }

  return totRate;
}



