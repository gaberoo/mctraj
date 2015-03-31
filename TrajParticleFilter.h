#ifndef __TRAJPARTICLEFILTER_H__
#define __TRAJPARTICLEFILTER_H__

#include "EpiState.h"
#include "TrajParticle.h"
#include "Model.h"

#include <omp.h>

#include <set>
using namespace std;

namespace MCTraj {
  typedef double (*TreeProbFun)(TrajParticle* traj, const vector<int>& nlin, 
                                double eventTime, void* pars);

  class TrajParticleFilter {
    public:
      TrajParticleFilter(const Model* m) 
        : model(m), curStep(0), curLin(1), tree(NULL), 
          logw(0.0), logp(0.0), vflag(0), simType(0)
      {
        pf.push_back(vector<TrajParticle>());
      }

      TrajParticleFilter(const TrajParticleFilter& x) 
        : model(x.model), curStep(x.curStep), 
          curLin(x.curLin), 
          pf(x.pf), tree(x.tree),
          etypeFun(x.etypeFun), logw(x.logw), logp(x.logp),
          vflag(x.vflag), simType(x.simType)
      {}

      virtual ~TrajParticleFilter() {}

      void push_back(const TrajParticle& tp) { pf[curStep].push_back(tp); }
      void setTree(const Tree* tree, int skip = 0, gsl_rng** rng = NULL);

      size_t stepTree(const void* pars, gsl_rng** rng, double dt = INFINITY);
      int stepAddTP(size_t j, const void* pars, gsl_rng* rng);
      size_t stepAdd(const void* pars, gsl_rng** rng);

      void calcWeights(const void* pars);
      int addTreeEvent(const void* pars, gsl_rng** rng, int noProb = 0);
      void setLast();

      int filter(gsl_rng** rng, char filter_type = 'w');
      int sampleInPlace(gsl_rng* rng);

      double meanWeight() const;
      void resetWeights();
      double meanCWeight() const;
      double w_vec(vector<double>& w) const;

      double meanProb() const;
      double p_vec(vector<double>& p) const;

      TrajParticle sample(gsl_rng* rng) const;

      size_t size() const { return pf[curStep].size(); }
      size_t curGen() const { return pf.size()-1; }

      TrajParticle& operator[](size_t i) { return pf[curStep][i]; }
      const Particle& at(size_t i) const { return pf[curStep].at(i); }

      void inc();

      double maxTime() const { return tree->maxTime(); }
      double tree_time() const { return tree->ttype(curStep); }
      int tree_ttype() const { return tree->ttype(curStep); }
      void printNames(char sep = '\n') const;

      void print(size_t j) const;
      void printAll() const;

      void set_etype_fun(void (*f)(int eventType, vector<int>& x)) { etypeFun= f; }

      double est() const { return logw; }
      double est_p() const { return logp; }

      void printFromLast(const void* pars = NULL, int offset = 0) const;
      void printFromFirst(const void* pars = NULL) const;

      void meanState(double mean[], size_t gen) const;
      void varState(double var[], const double mean[], size_t gen) const;

      void printMeanTraj(ostream& out) const;

      Trajectory singleTraj(gsl_rng* rng) const;

      size_t getCurStep() const { return curStep; }
      size_t getCurLin() const { return curLin; }

      void setVerbosity(int v) { vflag = v; }

    protected:
      const Model* model;

      inline TrajParticle& particle(size_t i) { return pf[curStep][i]; }
      inline const TrajParticle& particle(size_t i) const { return pf[curStep][i]; }

      size_t curStep;
      size_t curLin;
      vector< vector<TrajParticle> > pf;
      const Tree* tree;
      void (*etypeFun)(int eventType, vector<int>& x);
      double logw;
      double logp;
      int vflag;

      int simType;
  };
}

#endif
