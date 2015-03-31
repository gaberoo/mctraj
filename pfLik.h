#include "../Tree.h"

#include "TrajParticleFilter.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"
using namespace MCTraj;

#include <iostream>
#include <string>
#include <cstring>
using namespace std;

#include <gsl/gsl_sf.h>

namespace MCTraj {
  double pfLik(const Model* m, const EpiState& init, const Tree& tree,
               size_t num_particles, gsl_rng** rng, 
               double filter_time = 0.0,
               int vflag = 0, Trajectory* out = NULL, 
               int skip = 1, int print_particles = 0);
}

