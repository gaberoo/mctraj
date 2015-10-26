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

  typedef double (*RateFun)(const EpiState&, const void* pars, double& trueRate);
  typedef double (*ProbFun)(const EpiState&, const void* pars);
  typedef int (*BranchFun)(const EpiState& es, rng::RngStream* rng, 
                           StateTransition& st, const void* pars);

  typedef struct {
    size_t num_particles = 10;
    int print_particles = 0;
    int vflag = 0;
    int skip = 0;
    int reps = 1;
    int print_traj = 0;
    int full_tree = 0;
    int history = 0;
    char model_type = 'S';
    long unsigned seed = time(NULL);
    double filter_time = 0.0;
    double step_size = INFINITY;
    bool adj_zero = true;
  } PFPars;

  void pf_pars_read_json(PFPars* p, rapidjson::Value& json);
}


#endif // __MCTRAJ__
