#ifndef __HISTORYFILTER_H__
#define __HISTORYFILTER_H__

#include "TrajParticleFilter.h"

namespace MCTraj {
  class HistoryFilter : public TrajParticleFilter {
    public:
      HistoryFilter(const Model* m) 
        : TrajParticleFilter(m) 
      {
        pf.push_back(vector<TrajParticle>());
      }

      HistoryFilter(const HistoryFilter& x) 
        : TrajParticleFilter(x), pf(x.pf)
      {}

      virtual ~HistoryFilter() {}

      inline void push_back(const TrajParticle& tp) { pop().push_back(tp); }
      void setTree(const Tree* tree, int skip = 0, rng::Rng* rng = NULL);

      void setLast();
      void inc();

      void printMeanTraj(ostream& out) const;
      Trajectory singleTraj(rng::RngStream* rng) const;

      void printFromLast(const void* pars = NULL, int offset = 0) const;
      void printFromFirst(const void* pars = NULL) const;

      void meanState(double mean[], size_t gen) const;
      void varState(double var[], const double mean[], size_t gen) const;

      inline void reserve(size_t n) { pf.reserve(n); }
      inline size_t size() const { return pop().size(); }
      inline size_t curGen() const { return pf.size()-1; }

      void weights(vector<double>& mu, vector<double>& s2) const;

    protected:
      inline vector<TrajParticle>& pop() { return pf[curStep]; }
      inline vector<TrajParticle>& prev() { return pf[curStep-1]; }

      inline const vector<TrajParticle>& pop() const { return pf[curStep]; }
      inline const vector<TrajParticle>& prev() const { return pf[curStep-1]; }

      inline TrajParticle& particle(size_t i) { return pop()[i]; }
      inline const TrajParticle& particle(size_t i) const { return pop()[i]; }

      vector< vector<TrajParticle> > pf;
  };
}

#endif
