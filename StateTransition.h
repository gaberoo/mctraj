#ifndef __STATETRANSITION_H__
#define __STATETRANSITION_H__

#include "EpiState.h"
#include "BranchState.h"
#include "Trajectory.h"

namespace MCTraj {
  class StateTransition {
    friend class Trajectory;
    friend class EpiState;

    public:
      StateTransition() {}
      StateTransition(size_t nstates) 
        : time(0.0), prob(0.0), type(NULL), trans(nstates), 
          eventType(-1), absTime(0.0), relprob(0.0)
      {}

      StateTransition(const EpiState& es) 
        : time(0.0), prob(0.0), type(NULL), trans(es.numStates()), 
          eventType(-1), absTime(0.0), relprob(0.0)
      {}

      StateTransition(const StateTransition& st) 
        : time(st.time), prob(st.prob), type(st.type), 
          trans(st.trans),
          eventType(st.eventType), absTime(st.absTime), 
          branchTrans(st.branchTrans), 
          relprob(st.relprob)
      {}

      StateTransition(double t, const TransitionType& tt, const EpiState& es, 
                      const void* pars, int et = -1, double at = -INFINITY) 
        : time(t), prob(0.0), type(&tt), trans(tt.change), 
          eventType(et), absTime(at), relprob(0.0)
      { /* prob = tt.applyRate(es,pars); */ }

      virtual ~StateTransition() {}

      int& operator[](size_t i) { return trans[i]; }
      int operator[](size_t i) const { return trans[i]; }

      double atTime() const { return time; }
      double realTime() const { return absTime; }

      double getProb() const { return prob; }

      const TransitionType* getType() const { return type; }

      const vector<int>& getTrans() const { return trans; }

      int alters(size_t i) const { return type->alters(i); }

      int etype() const { return eventType; }
      void etype(int et) { eventType = et; }

      StateTransition& operator=(const StateTransition& st) {
        time = st.time;
        prob = st.prob;
        type = st.type;
        trans = st.trans;
        eventType = st.eventType;
        absTime = st.absTime;
        branchTrans = st.branchTrans;
        relprob = st.relprob;
        return *this;
      }

      inline void addBranchTrans(int id, int old_color, int new_color) {
        branchTrans.push_back(BranchStateChange(id,old_color,new_color));
      }

      friend ostream& operator<<(ostream& out, const StateTransition& trans);

      string to_json() const;

    protected:
      double time;
      double prob;
      const TransitionType* type;

      vector<int> trans;

      int eventType;
      double absTime;

    public:
      vector<BranchStateChange> branchTrans;
      double relprob;
  };
}

#endif
