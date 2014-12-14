#ifndef __PARTICLEFILTER_H__
#define __PARTICLEFILTER_H__

#include <vector>
#include <string>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Particle.h"

namespace MCTraj {
  class ParticleFilter {
    public:
      ParticleFilter() {}
      ParticleFilter(const ParticleFilter& pf) : pf(pf.pf) {}
      virtual ~ParticleFilter() {}

      void filter(gsl_rng* rng);
      const Particle& at(size_t i) const { return pf.at(i); }

    protected:
      vector<Particle> pf;
  };
}

#endif
