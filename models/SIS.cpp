#include "SIS.h"

double MCTraj::SISModel::infRateFun(const EpiState& es, const void* pars) 
{
  // cerr << "SIS:INF" << endl;
  EpiPars ep = *(EpiPars*) pars;
  double lambdaI = ep.beta*es[0]*es[1]/ep.N;
  return (lambdaI > 0.0) ? lambdaI : 0.0;
}

double MCTraj::SISModel::treeProbInf(const EpiState& es, const void* pars) 
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

double MCTraj::SISModel::recovRateFun(const EpiState& es, const void* pars) 
{
  // cerr << "SIS:RECOV" << endl;
  EpiPars ep = *(EpiPars*) pars;
  double I = es[1];
  return (ep.mu+ep.psi)*I;
}

double MCTraj::SISModel::treeProbRecov(const EpiState& es, const void* pars) 
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

double MCTraj::SISModel::treeObsInf(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double S = es[0];
  double I = es[1];
  double lambda = ep.beta*S/ep.N;
  return 2.0*lambda/(I+1.);
  // return lambda;
}

/************************************************************************/

double MCTraj::SISModel::treeObsRecov(const EpiState& es, const void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[1];
  return (ep.psi > 0) ? ep.psi*I : 1.0;
}

/************************************************************************/

double MCTraj::SIS::sample_rho(const EpiState& es, gsl_rng* rng, void* pars) const {
  int k = es[2];
  int I = es[1];
  double w = 0.0;
  if (rho >= 1.0 && I == k) {
    w = 1.0;
  } else if (rho < 1.0 && I >= k) {
    w = gsl_ran_binomial_pdf(k,rho,I);
  }
  return w;
}

