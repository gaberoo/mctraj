#include "SIR.h"

double MCTraj::SIRModel::infRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double lambdaI = ep.beta*es[0]*es[1]/ep.N;
  return (lambdaI > 0.0) ? lambdaI : 0.0;
}

double MCTraj::SIRModel::treeProbInf(const EpiState& es, const void* pars) 
{
  /* Three possibilities upon infection:
   *   (i)   infection involving two tree branches
   *   (ii)  infection involving one tree branch
   *   (iii) infeciton involving zero tree branches
   * The probability that it doesn't involve two tree branches is,
   *   p = 1 - (k choose 2) / (I choose 2) 
   */
  double I = es[1];
  double k = es[2];
  double p = 1.0 - k*(k-1.)/(I*(I-1.));
  return (I >= k) ? p : 0.0;
}

/************************************************************************/

double MCTraj::SIRModel::recovRateFun(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[1];
  return (ep.mu+ep.psi)*I;
}

double MCTraj::SIRModel::treeProbRecov(const EpiState& es, const void* pars) 
{
  /* Four possibilities:
   *   (i)   recovery in a non-tree branch
   *   (ii)  sampling in a non-tree branch
   *   (iii) recovery in a tree branch
   *   (iv)  sampling in a tree branch 
   * The probability that it's a recovery in a non-tree branch is
   *   p = (1-s) 
   */
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.mu);
  return (es[1]+1 > es[2]) ? (1.-s) : 0.0;
}

/************************************************************************/

double MCTraj::SIRModel::treeObsInf(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double S = es[0];
  double I = es[1];
  double lambda = ep.beta*S/ep.N;
  return 2.0*lambda/(I+1.);
  // return lambda;
}

/************************************************************************/

double MCTraj::SIRModel::treeObsRecov(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[1];
  return (ep.psi > 0) ? ep.psi*I : 1.0;
}

