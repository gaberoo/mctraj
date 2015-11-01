#ifndef __TRAJPARTICLEFILTER_H__
#define __TRAJPARTICLEFILTER_H__

#include <rng/Rng.h>
#include "EpiState.h"
#include "TrajParticle.h"
#include "Model.h"
#include "ascii.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <set>
using namespace std;

#include <gsl/gsl_statistics.h>

namespace MCTraj {
  typedef double (*TreeProbFun)(TrajParticle* traj, const vector<int>& nlin, 
                                double eventTime, void* pars);

  class TrajParticleFilter {
    public:
      TrajParticleFilter(const Model* m) 
        : model(m), curStep(0), curLin(1),
          tree(NULL), logw(0.0), logp(0.0), vflag(0), simType(0)
      {}

      TrajParticleFilter(const TrajParticleFilter& x) 
        : model(x.model), curStep(x.curStep), 
          curLin(x.curLin), tree(x.tree),
          etypeFun(x.etypeFun), logw(x.logw), logp(x.logp),
          vflag(x.vflag), simType(x.simType)
      {}

      virtual ~TrajParticleFilter() {}

      virtual void push_back(const TrajParticle& tp) = 0;
      virtual void setTree(const Tree* tree, int skip = 0, rng::Rng* rng = NULL);

      size_t stepTree(const Pars* pars, rng::Rng* rng, bool adjZero = true, double dt = INFINITY);
      void calcWeights(const void* pars);
      int addTreeEvent(const void* pars, rng::Rng* rng);

      virtual void setLast() = 0;
      virtual void inc() = 0;

      int filter(rng::Rng* rng, char filter_type = 'w');
      int sampleInPlace(rng::RngStream* rng);

      double meanWeight() const;
      double meanProb() const;
      double w_vec(vector<double>& w) const;
      double p_vec(vector<double>& p) const;

      TrajParticle sample(rng::RngStream* rng) const;

      virtual size_t size() const = 0;
      virtual size_t curGen() const = 0;

      TrajParticle& operator[](size_t i) { return particle(i); }
      const Particle& at(size_t i) const { return particle(i); }

      inline double cur_time() const { return particle(0).getTime(); }

      // tree properties
      double maxTime() const { return tree->maxTime(); }
      double tree_time() const { return tree->ttype(curStep); }
      int tree_ttype() const { return tree->ttype(curStep); }

      // print trajectories
      void printNames(char sep = '\n') const;
      void print(size_t j) const;
      void printAll() const;

      inline void set_etype_fun(void (*f)(int eventType, vector<int>& x)) { etypeFun= f; }
      inline void setVerbosity(int v) { vflag = v; }

      inline double est() const { return logw; }
      inline double est_p() const { return logp; }

      virtual void meanState(double mean[], size_t gen) const = 0;
      virtual void varState(double var[], const double mean[], size_t gen) const = 0;

      inline size_t getCurStep() const { return curStep; }
      inline size_t getCurLin() const { return curLin; }

      inline const TransitionType* getCurType() const {
        return model->getObsType(model->mapType(tree_ttype()));
      }

      virtual void copyFromPrev(size_t i, size_t j) = 0;

      // extract trajectories
      virtual Trajectory singleTraj(int traj_id) const = 0;

      inline const TrajParticle& cparticle(size_t i) const { return particle(i); }

    protected:
      const Model* model;

      virtual TrajParticle& particle(size_t i) = 0;
      virtual const TrajParticle& particle(size_t i) const = 0;

      virtual vector<TrajParticle>& pop() = 0;
      virtual const vector<TrajParticle>& pop() const = 0;
      virtual vector<TrajParticle>& prev() = 0;
      virtual const vector<TrajParticle>& prev() const = 0;

      const TrajParticle& prev(size_t i) const { return prev()[i]; }

      size_t curStep;
      size_t curLin;

      const Tree* tree;
      void (*etypeFun)(int eventType, vector<int>& x);
      double logw;
      double logp;
      int vflag;

      int simType;

      const PFPars* pars;
  };
}

#endif
