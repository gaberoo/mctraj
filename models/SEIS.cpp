#include "SEIS.h"

double MCTraj::SEISModel::infRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double lambdaI = ep.beta*es[0]*es[2]/ep.N;
  // double lambdaI = ep.beta*es[2];
  return (lambdaI > 0.0) ? lambdaI : 0.0;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::infTreeProb(const EpiState& es, const void* pars) 
{
  /* Three possibilities upon infection:
   *   (i)   infection involving two tree branches
   *   (ii)  infection involving one tree branch
   *   (iii) infeciton involving zero tree branches
   * The probability that it doesn't involve two tree branches is,
   *   p = 1 - (k choose 2) / (I choose 2) 
   */
  double E = es[1];
  double I = es[2];
  double kE = es[3];
  double kI = es[4];

  // double p = 1.0 - kI/I * kE/E;
  double p = 1.0;

#ifdef DEBUG
  cerr << "treeProbInf   :: ES = " << es << " | p = " << p << endl;
#endif
  return (I >= kI && E >= kE) ? p : 0.0;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::infValidate(const EpiState& es, const void* pars) {
  return 1;
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::infTreeObs(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double S = es[0];
  // double E = es[1];
  double I = es[2];
  double lambda = ep.beta*S/ep.N;
  // double lambda = ep.beta;
  int branch_color = 0;
  double x = 0.0;
  if (I <= 0.0) {
#ifdef DEBUG
    cerr << "treeObsInf    :: ES = " << es << " | no infecteds!" << endl;
#endif
    return 0.0;
  }
  if (es.curBranch.size() > 0) {
    branch_color = es.branches.getCol(es.curBranch[0]);
//    x = (branch_color == 1) ? lambda/(E+1.) : 0.0; // prob of the particular branch
    x = (branch_color == 1) ? lambda : 0.0; // prob of the particular branch
#ifdef DEBUG
    if (x <= 0.0) {
      cerr << "treeObsInf    :: ES = " << es 
           << " | " << es.curBranch[0] << " -- " << branch_color << endl;
    }
#endif
    return x;
  } else {
#ifdef DEBUG
    cerr << "treeObsInf    :: no node information! ES = " << es << endl;
#endif
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

  double p = 1.0;
  double dt = INFINITY;

  // printf(">>>  I = %.0f, kI = %.0f\n",I,kI);
  // the probability of picking a tree branch is kI/I
  if (r < kI/I) {
    // get a random branch in state I (I=1)
    int id = es.branches.random_color(rng,1);
    if (id >= 0) {
      if (ep->tree != NULL) {
//        cerr << id << "/" << ep->tree->branches.size() << endl;
//        cerr << "Infection happened on branch "
//            << ep->tree->branches[id].id << ". " << endl;
//        cerr << "   Time left before next event = "
//            << ep->tree->branches[id].time - es.time << endl;
//        cerr << "   Type of next event = "
//            << ep->tree->branches[id].type << endl;

        // the probability that the new individual becomes this 
        // branch is 0.5 * (1 + exp(-dt)), where dt is the time 
        // to the next interval
        dt = ep->tree->branches[id].time - es.time;
        switch (ep->tree->branches[id].type) {
          case 0:
          case 1:
          case 2: p = 1.0 + exp(-dt); break;
          default: p = 1.0; break;
        }
      }

      rng->uniform(1,&r);
      if (r < 0.5*p) {
        // the color change happens on the branch
        st.branchTrans.push_back(BranchStateChange(id,1,0));
        st[3] = 1;
        st[4] = -1;
        st.relprob = 1.0/p;
      } else {
        st.relprob = 1.0/(2.0-p);
      }
    } else {
      // otherwise raise an error
      cerr << "[I] No branch of color 'I' exists! | ES = " << es << endl;
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

  int b1 = es.curBranch.at(1);
  int b2 = es.curBranch.at(2);

  double dt1 = ep->tree->branches[b1].time - es.time;
  double dt2 = ep->tree->branches[b2].time - es.time;

  double p = dt1/(dt1+dt2);

  int a, b;
  double r;
  rng->uniform(1,&r);
  if (r < p) { 
    a = 1; 
    b = 2; 
    st.relprob = 0.5/p;
  } else { 
    b = 1; 
    a = 2; 
    st.relprob = 0.5/(1.0-p);
  }

  // get colors of the branches
  int acol = es.branches.getCol(es.curBranch.at(a));
  int bcol = es.branches.getCol(es.curBranch.at(b));

//  cerr << "infection observed at ("
//       << es.curBranch.at(0) << "," 
//       << es.curBranch.at(1) << ","
//       << es.curBranch.at(2) << ")" << endl;

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
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(a),acol,1));
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(b),bcol,0));
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
    }
  }
  return 0;
}

/****************************************************************************/
/* RECOVERY FUNCTIONS *******************************************************/
/****************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  return (ep.mu+ep.psi)*es[2];
}

//----------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeProb(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.mu);
  double I = es[2]+1.0;
  double kI = es[4];
  double dw = (I > kI) ? (1.-s) : 0.0;
  if (es[2] < es[4] || es[1] < es[3]) dw = 0.0;
#ifdef DEBUG
  cerr << "treeProbRecov :: ES = " << es << " | w = " << dw 
       << " >> " << I << " --> " << kI << endl;
#endif
  return dw;
}

//----------------------------------------------------------------------------

int MCTraj::SEISModel::recovValidate(const EpiState& es, const void* pars) {
  return 1;
}

//---------------------------------------------------------------------------

double MCTraj::SEISModel::recovTreeObs(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  // double I = es[2];
  int branch_color = 0;
  if (es.curBranch.size() > 0) {
    branch_color = es.branches.getCol(es.curBranch[0]);
//    double x = (branch_color == 1) ? ep.psi*I : 0.0;
    double x = (branch_color == 1) ? ep.psi : 0.0;
#ifdef DEBUG
    if (x <= 0.0) {
      cerr << "treeObsRecov :: ES = " << es 
           << " | " << es.curBranch[0] << " -- " << branch_color << endl;
    }
#endif
    return x;
  } else {
    cerr << "treeObsRecov :: no node information! ES = " << es << endl;
    return 1.0;
  }
}

/************************************************************************/

int MCTraj::SEISModel::recovBranchObs(const EpiState& es, rng::RngStream* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
//    cerr << ">> " << es.curBranch.at(0) << endl;
    st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
  }
  return 0;
}

/* TRANSITIONS **********************************************************/

double MCTraj::SEISModel::transRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double E = es[1];
  return ep.gamma*E;
}

double MCTraj::SEISModel::transTreeProb(const EpiState& es, const void* pars) 
{
  double E = es[1];
  double kE = es[3];
  // return (E >= kE) ? (E-kE)/E : 0.0;
  return (E >= kE) ? 1.0 : 0.0;
}

int MCTraj::SEISModel::validateTrans(const EpiState& es, const void* pars) {
  return 1;
}

double MCTraj::SEISModel::transTreeObs(const EpiState& es, const void* pars) 
{
  debug("treeObsTrans :: ES = (%d,%d,%d,%d,%d)",es[0],es[1],es[2],es[3],es[4]);
  return 1.0;
}

/************************************************************************/

int MCTraj::SEISModel::transBranch(const EpiState& es, rng::RngStream* rng, 
                                   StateTransition& st, const void* pars)
{
  double E = es[1];
  double kE = es[3];
  double r;
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
  return 0;
}

/************************************************************************/

int MCTraj::SEISModel::transBranchObs(const EpiState& es, rng::RngStream* rng, 
                                      StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
#ifdef DEBUG
    cerr << ">> " 
         << es.curBranch.at(0) << ", "
         << es.curBranch.at(1) << ", "
         << es.curBranch.at(2) << " () " 
         << es.branches.getCol(es.curBranch.at(0)) << endl;
#endif
    if (es.curBranch.at(0) >= 0) {
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),0,-1));
    }
    if (es.curBranch.at(1) >= 0) {
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(1),-1,1));
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


