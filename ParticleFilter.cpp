#include "ParticleFilter.h"

namespace MCTraj {
  void ParticleFilter::filter(gsl_rng* rng) {
    vector<double> w(pf.size(),0.0);
    vector<Particle> old(pf);
    for (size_t i(0); i < pf.size(); ++i) w[i] = pf[i].getWeight();
    gsl_ran_discrete_t* pp = gsl_ran_discrete_preproc(pf.size(),w.data());
    for (size_t i(0); i < pf.size(); ++i) {
      size_t j = gsl_ran_discrete(rng,pp);
      pf[i] = old.at(j);
    }
    gsl_ran_discrete_free(pp);
  }
}
