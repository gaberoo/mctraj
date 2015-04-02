#ifndef __MCTRAJ__
#define __MCTRAJ__

#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

#include <rng/Rng.h>
#include "ParticleFilter.h"

#ifndef DEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>

namespace MCTraj {
  class EpiState;
  class TransitionType;
  class StateTransition;
  class Trajectory;
  class Model;

  typedef double (*RateFun)(const EpiState&, const void* pars);
  typedef double (*ProbFun)(const EpiState&, const void* pars);
  typedef int (*BranchFun)(const EpiState& es, rng::RngStream* rng, 
                           StateTransition& st, const void* pars);

  typedef struct {
    size_t num_particles;
    int print_particles;
    int vflag;
    int skip;
    int reps;
    int model_type;
    int print_traj;
    int full_tree;
    long unsigned seed;
    double filter_time;
    bool adj_zero;
  } PFPars;
}

#endif // __MCTRAJ__
