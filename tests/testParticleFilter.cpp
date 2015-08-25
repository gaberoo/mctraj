#include <iostream>
#include <string>
#include <cstring>
using namespace std;

#include "ParticleFilter.h"
using namespace MCTraj;

int main() {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,time(NULL));

  ParticleFilter pf;
  for (size_t i(0); i < 10; ++i) {
    char name[32];
    sprintf(name,"P%d",(int) i);
    pf.push_back(Particle(string(name),gsl_rng_uniform(rng)));
  }

  for (size_t i(0); i < pf.size(); ++i) {
    cout << i << " " << pf[i].getName() << " " << pf[i].getWeight() << endl;
  }

  pf.filter(rng);

  for (size_t i(0); i < pf.size(); ++i) {
    cout << i << " " << pf[i].getName() << " " << pf[i].getWeight() << endl;
  }

  gsl_rng_free(rng);
  return 0;
}

