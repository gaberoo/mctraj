#ifndef __MCTRAJ__
#define __MCTRAJ__

#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "ParticleFilter.h"

namespace MCTraj {
  class EpiState;
  class TransitionType;
  class StateTransition;
  class Trajectory;
  class Model;

  typedef double (*RateFun)(const EpiState&, const void* pars);
  typedef double (*ProbFun)(const EpiState&, const void* pars);
  typedef void (*BranchFun)(const EpiState& es, gsl_rng* rng, 
                            StateTransition& st, const void* pars);
}

#endif // __MCTRAJ__
