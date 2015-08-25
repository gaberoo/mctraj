#include <iostream>
#include <ctime>
using namespace std;

#include <rng/GSLStream.h>

#include "MCTraj.h"
#include "EpiState.h"
#include "Trajectory.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"
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
  rng::GSLStream rng;
  rng.alloc(1337);

  SEISModel::EpiPars seis_pars;
  seis_pars.N = 100;
  seis_pars.beta = 1.0;
  seis_pars.mu = 0.1;
  seis_pars.psi = 0.1;
  seis_pars.rho = 0.0;
  seis_pars.gamma = 0.01;
  Model* mpt = new SEIS(&seis_pars);

  EpiState* es = new EpiState(SEISModel::nstates);
  (*es)[0] = ((int) seis_pars.N)-6;
  (*es)[1] = 3;
  (*es)[2] = 3;
  (*es)[3] = 1;
  (*es)[4] = 2;

  vector<double> transRates(mpt->ntrans(),0.0);
  mpt->calculateTransRates(*es,transRates);


  double r = 0.0;
  for (size_t i = 0; i < transRates.size(); ++i) {
    r = mpt->getTType(i)->applyRate(*es,&seis_pars);
    cout << i << " " << r/transRates.back() << " " << transRates[i] << endl;
  }

  cout << endl;

  size_t ti; 
  size_t ntot = 1000000;
  vector<int> cnt(mpt->ntrans());
  for (size_t i = 0; i < ntot; ++i) {
    ti = mpt->chooseTransition(&rng,transRates);
    cnt[ti]++;
  }
  for (size_t i = 0; i < cnt.size(); ++i) {
    cout << i << " " << cnt[i] << " " << 1.*cnt[i]/ntot << endl;
  }

//  Trajectory T(es);
//  T.addTransType(nstates,infChange,infRateFun);
//  T.addTransType(nstates,recoverChange,recovRateFun);
//
//  T.simulateTrajectory(30.0,&pars,rng);
//  cout << T << endl;

  return 0;
}

