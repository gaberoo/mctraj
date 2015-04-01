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
               const PFPars& pars, Rng* rng, Trajectory* out = NULL);
}

