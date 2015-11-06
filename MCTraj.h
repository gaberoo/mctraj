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

  class PFPars {
    public:
      PFPars()
        : num_particles(10), print_particles(0),
          vflag(0), skip(0), reps(1), print_traj(0),
          full_tree(0), history(0), model_type('S'),
          seed(time(NULL)), filter_time(0.0), 
          step_size(INFINITY), adj_zero(true)
      {}
      virtual ~PFPars() {}

      void from_json(rapidjson::Value& json);

    public:
      size_t num_particles;
      int print_particles;
      int vflag;
      int skip;
      int reps;
      int print_traj;
      int full_tree;
      int history;
      char model_type;
      long unsigned seed;
      double filter_time;
      double step_size;
      bool adj_zero;
  };

  inline void pf_pars_read_json(PFPars* p, rapidjson::Value& json) {
    p->from_json(json);
  }
}


#endif // __MCTRAJ__
