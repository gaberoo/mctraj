#include "SEIS.h"

double MCTraj::SEISModel::infRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double lambdaI = ep.beta*es[0]*es[2]/ep.N;
  return (lambdaI > 0.0) ? lambdaI : 0.0;
}

/************************************************************************/

double MCTraj::SEISModel::treeProbInf(const EpiState& es, const void* pars) 
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

  // double p = 1.0 - k*(k-1.)/(I*(I-1.));
  double p = 1.0 - kI*kE/(E*I);

#ifdef DEBUG
  cerr << "treeProbInf   :: ES = " << es << " | p = " << p << endl;
#endif
  return (I >= kI && E >= kE) ? p : 0.0;
}

/************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  return (ep.mu+ep.psi)*es[2];
}

/************************************************************************/

double MCTraj::SEISModel::treeProbRecov(const EpiState& es, const void* pars) 
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

/************************************************************************/

double MCTraj::SEISModel::treeObsInf(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double S = es[0];
  double E = es[1];
  double I = es[2];
  double lambda = ep.beta*S/ep.N;
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
    x = (branch_color == 1) ? lambda/(E+1.) : 0.0; // prob of the particular branch
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

/************************************************************************/

double MCTraj::SEISModel::treeObsRecov(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[2];
  int branch_color = 0;
  if (es.curBranch.size() > 0) {
    branch_color = es.branches.getCol(es.curBranch[0]);
    double x = (branch_color == 1) ? ep.psi*I : 0.0;
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

double MCTraj::SEISModel::transRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double E = es[1];
  return ep.gamma*E;
}

double MCTraj::SEISModel::treeProbTrans(const EpiState& es, const void* pars) 
{
  double E = es[1];
  double kE = es[3];
  // return (E >= kE) ? (E-kE)/E : 0.0;
  return (E >= kE) ? 1.0 : 0.0;
}

double MCTraj::SEISModel::treeObsTrans(const EpiState& es, const void* pars) 
{
  debug("treeObsTrans :: ES = (%d,%d,%d,%d,%d)",es[0],es[1],es[2],es[3],es[4]);
  return 1.0;
}

/************************************************************************/

int MCTraj::SEISModel::branchInf(const EpiState& es, rng::RngStream* rng, 
                                 StateTransition& st, const void* pars)
{
  // double E = es[1];
  double I = es[2];
  // double kE = es[3];
  double kI = es[4];
  double r; 
  rng->uniform(1,&r);
  // printf(">>>  I = %.0f, kI = %.0f\n",I,kI);
  if (r < 0.5*kI/I) {
    int id = es.branches.random_color(rng,1);
    if (id >= 0) {
      st.branchTrans.push_back(BranchStateChange(id,1,0));
      st[3] = 1;
      st[4] = -1;
    } else {
      // cerr << "[I] No branch of color 'I' exists! | ES = " << es << endl;
      // cerr << es.branches.to_json() << endl;
      return -1;
    }
  }
  return 0;
}

/************************************************************************/

int MCTraj::SEISModel::branchTrans(const EpiState& es, rng::RngStream* rng, 
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

int MCTraj::SEISModel::obsBranchInf(const EpiState& es, rng::RngStream* rng, 
                                     StateTransition& st, const void* pars)
{
  int r;
  rng->uniform_int(1,&r,0,2);
  int a, b;
  if (r) { a = 1; b = 2; } else { b = 1; a = 2; }
  int acol = es.branches.getCol(es.curBranch.at(a));
  int bcol = es.branches.getCol(es.curBranch.at(b));
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

/************************************************************************/

int MCTraj::SEISModel::obsBranchRecov(const EpiState& es, rng::RngStream* rng, 
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

/************************************************************************/

int MCTraj::SEISModel::obsBranchTrans(const EpiState& es, rng::RngStream* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
//    cerr << ">> " 
//         << es.curBranch.at(0) << ", "
//         << es.curBranch.at(1) << ", "
//         << es.curBranch.at(2) << " () " 
//         << es.branches.getCol(es.curBranch.at(0)) << endl;
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

