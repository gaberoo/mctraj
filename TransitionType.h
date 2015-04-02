#ifndef __TRANSITIONTYPE_H__
#define __TRANSITIONTYPE_H__

#include <rng/Rng.h>
#include "MCTraj.h"
#include "EpiState.h"

namespace MCTraj {
  class TransitionType {
    friend class StateTransition;
    friend class EpiState;

    public:
      TransitionType() : name("") {}

      TransitionType(const TransitionType& tt) 
        : name(tt.name), change(tt.change), rate(tt.rate), prob(tt.prob)
      {}

      TransitionType(size_t nstates, const int* changeArray, 
                     RateFun r = NULL, ProbFun p = NULL, BranchFun b = NULL) 
        : name(""), change(nstates), rate(r), prob(p), bfun(b)
      { change.assign(changeArray,changeArray+nstates); }

      TransitionType(string name, size_t nstates, const int* changeArray, 
                     RateFun r = NULL, ProbFun p = NULL, BranchFun b = NULL) 
        : name(name), change(nstates), rate(r), prob(p), bfun(b)
      { change.assign(changeArray,changeArray+nstates); }

      virtual ~TransitionType() {}

      inline double applyRate(const EpiState& es, const void* pars) const 
      { if (rate != NULL) { return rate(es,pars); } else { return 0.0; } }

      inline double applyProb(const EpiState& es, const void* pars) const 
      { if (prob != NULL && pars != NULL) { return prob(es,pars); } else { return 0.0; } }

      inline int applyBranch(const EpiState& es, rng::RngStream* rng, 
                             StateTransition& st, const void* pars) const 
      { return (bfun != NULL) ? bfun(es,rng,st,pars) : -100; }

      inline int alters(size_t i) const { return change.at(i); }
      inline int operator[](size_t i) const { return change[i]; }

      friend ostream& operator<<(ostream& out, const TransitionType& tt);

      inline TransitionType& operator=(const TransitionType& tt) {
        change = tt.change;
        rate = tt.rate;
        prob = tt.prob;
        return *this;
      }

      inline void setName(string n) { name = n; }
      inline string getName() const { return name; }

    protected:
      string name;
      vector<int> change;
      RateFun rate;
      ProbFun prob;
      BranchFun bfun;
  };
}

#endif

