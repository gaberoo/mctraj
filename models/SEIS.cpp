#include "SEIS.h"

double MCTraj::SEISModel::infRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double lambdaI = ep.beta*es[0]*es[2]/ep.N;
  return (lambdaI > 0.0) ? lambdaI : 0.0;
}

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

  return (I >= kI && E >= kE) ? p : 0.0;
}

/************************************************************************/

double MCTraj::SEISModel::recovRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[2];
  return (ep.mu+ep.psi)*I;
}

double MCTraj::SEISModel::treeProbRecov(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.mu);
  double I = es[2]+1.0;
  double kI = es[4];
  return (I > kI) ? (1.-s) : 0.0;
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
  if (I <= 0.0) return 0.0;
  if (es.curBranch.size() > 0) {
    branch_color = es.branches.getCol(es.curBranch[0]);
    x = (branch_color == 1) ? lambda/(E+1.) : 0.0;
//    cerr << "treeObsInf :: ES = " << es 
//         << " | " << es.curBranch[0] << " -- " << branch_color << endl;
    return x;
  } else {
    cerr << "treeObsInf :: no node information! ES = " << es << endl;
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
//    cerr << "treeObsRecov :: ES = " << es 
//         << " | " << es.curBranch[0] << " -- " << branch_color 
//         << "(( " << x << endl;
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
  return 1.0;
}

/************************************************************************/

void MCTraj::SEISModel::branchInf(const EpiState& es, gsl_rng* rng, 
                                  StateTransition& st, const void* pars)
{
  // double E = es[1];
  double I = es[2];
  // double kE = es[3];
  double kI = es[4];
  double r = gsl_rng_uniform(rng);
  if (r < 0.5*kI/I) {
    int id = es.branches.random_color(rng,1);
    if (id >= 0) {
      st.branchTrans.push_back(BranchStateChange(id,1,0));
      st[3] = 1;
      st[4] = -1;
    } else {
      cerr << "[I] Color doesn't exist (" << id << ") !" << endl;
    }
  }
}

void MCTraj::SEISModel::branchTrans(const EpiState& es, gsl_rng* rng, 
                                    StateTransition& st, const void* pars)
{
  double E = es[1];
  double kE = es[3];
  double r = gsl_rng_uniform(rng);
  if (r < kE/E) {
    int id = es.branches.random_color(rng,0);
    // cerr << "Trans on branch. " << r << "/" << id << " " << es << endl;
    if (id >= 0) {
      st.branchTrans.push_back(BranchStateChange(id,0,1));
      st[3] = -1;
      st[4] = 1;
    } else {
      cerr << "[T] Color doesn't exist (" << id << ") !" << endl;
    }
  }
}

/************************************************************************/

void MCTraj::SEISModel::obsBranchInf(const EpiState& es, gsl_rng* rng, 
                                     StateTransition& st, const void* pars)
{
  int r = gsl_rng_uniform_int(rng,2);
  int a, b;
//  int branch_color = es.branches.getCol(es.curBranch.at(0));
  if (r) { a = 1; b = 2; } else { b = 1; a = 2; }
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
//      cerr << es.curBranch.at(a) << " >> "
//           << es.curBranch.at(b) << " || "
//           << es.curBranch.at(0) << endl;
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(a),-1,1));
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(b),-1,0));
      st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
    }
  }
}

void MCTraj::SEISModel::obsBranchRecov(const EpiState& es, gsl_rng* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
//    cerr << ">> " << es.curBranch.at(0) << endl;
    st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),1,-1));
  }
}

void MCTraj::SEISModel::obsBranchTrans(const EpiState& es, gsl_rng* rng, 
                                       StateTransition& st, const void* pars)
{
  if (es.curBranch.size() < 1) {
    cerr << "Need node information to color tree!" << endl;
  } else {
//    cerr << ">> " << es.curBranch.at(0) << " () " 
//         << es.branches.getCol(es.curBranch.at(0)) << endl;
    st.branchTrans.push_back(BranchStateChange(es.curBranch.at(0),0,-1));
    st.branchTrans.push_back(BranchStateChange(es.curBranch.at(1),-1,1));
  }
}


