#include <iostream>
#include <ctime>
using namespace std;

#include "MCTraj.h"
using namespace MCTraj;

struct EpiPars {
  double beta;
  double gamma;
  double N;
};

double infRateFun(const EpiState& es, void* pars) {
  EpiPars ep = *(EpiPars*) pars;
  return ep.beta*es[0]*es[1]/ep.N;
}

double recovRateFun(const EpiState& es, void* pars) {
  EpiPars ep = *(EpiPars*) pars;
  return ep.gamma*es[1];
}

int main() {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,time(NULL));

  size_t nstates = 3;

  int infChange[] = { -1, 1, 0 };
  int recoverChange[] = { 0, -1, 1 };

  EpiPars pars = { 1.0, 0.1, 100.0 };

  EpiState es(nstates);
  es[0] = 99;
  es[1] = 1;

  Trajectory T(es);
  T.addTransType(nstates,infChange,infRateFun);
  T.addTransType(nstates,recoverChange,recovRateFun);

  T.simulateTrajectory(30.0,&pars,rng);
  cout << T << endl;

  gsl_rng_free(rng);
  return 0;
}

