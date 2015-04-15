#ifndef __STATICFILTER_H__
#define __STATICFILTER_H__

#include <iostream>
#include <vector>
using namespace std;

#include "TrajParticleFilter.h"

namespace MCTraj {
  typedef vector<TrajParticle> TPArr;

  class StaticFilter : public TrajParticleFilter {
    public:
      StaticFilter(const Model* m) 
        : TrajParticleFilter(m)
      { 
        A = &pf1; 
        B = &pf2; 
      }

      StaticFilter(const StaticFilter& x) : TrajParticleFilter(x) {}
      virtual ~StaticFilter();

      inline void push_back(const TrajParticle& tp) { 
        A->push_back(tp); 
        B->push_back(tp); 
      }
      inline void setLast() {}
      inline void inc() {
        ++curStep;
        if (B->size() != A->size()) B->resize(A->size());
        swap();
      }

      inline size_t size() const { return A->size(); }
      inline size_t curGen() const { return 0; }

      void meanState(double mean[], size_t gen = 0) const;
      void varState(double var[], const double mean[], size_t gen = 0) const;

      inline void copyAB() { *B = *A; }
      void copyFromPrev(size_t i, size_t j);

    protected:
      inline void swap() { TPArr* tmp = A; A = B; B = tmp; }

      inline TrajParticle& particle(size_t i) { return (*A)[i]; }
      inline const TrajParticle& particle(size_t i) const { return (*A)[i]; }

      inline vector<TrajParticle>& pop() { return *A; }
      inline vector<TrajParticle>& prev() { return *B; }

      inline const vector<TrajParticle>& pop() const { return *A; }
      inline const vector<TrajParticle>& prev() const { return *B; }

      TPArr pf1;
      TPArr pf2;

      TPArr* A;
      TPArr* B;
  };
}

#endif
