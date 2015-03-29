#include "TrajParticle.h"

namespace MCTraj {
  int TrajParticle::lastType(int type) const {
    if (nTrans() > 0) {
      for (int i(nTrans()-1); i > getLast(); --i) {
        if (transitions[i].etype() == type) {
          return i;
        }
      }
    }
    return -1;
  }

  // =========================================================================

  double TrajParticle::calcWeight(const void* pars) {
    double dw = trajProb(pars,getLast());
    weight *= dw;
    return weight;
  }

  // =========================================================================

  void TrajParticle::print_from_last(ostream& out, string prefix, 
                                     const void* pars, int offset) const {
    EpiState curState = getState();
    StateTransition st;
    const TransitionType* tt;

    double time = getTime();
    // int last = getLast();
    int last = getLast();
    if (offset > 0 && offset < (int) ivals.size()) {
      last = *(ivals.rbegin() + offset);
    }

    double w = 0.0;

    out << prefix << " "
        << name << " " 
        << id << " " << parent << " " 
        << setw(12) << getWeight() << " "
        << setw(12) << time << " " 
        << curState << " " << last << endl;

    if (transitionCount() > 0) {
      int i;
      for (i = transitionCount()-1; i > last; --i) {
        st = getTrans(i);
        tt = st.getType();
        w = tt->applyProb(curState,pars);
        time -= st.atTime();
        curState -= st.getTrans();
        out << prefix << " "
            << name << " " 
            << id << " " << parent << " " 
            << setw(12) << w << " "
            << setw(12) << time << " " 
            << curState << " " << i << endl;
      }
    }
  }

  // =========================================================================

  void TrajParticle::print_from_first(ostream& out, string prefix, const void* pars) const 
  {
    EpiState curState(initialState);
    StateTransition st;
    const TransitionType* tt = NULL;
    double w = 0.0;

    double time = init_time;

    out << prefix << " "
        << name << " " 
        << id << " " << parent << " " << getWeight() << " "
        << setw(12) << time << " " 
        << curState << endl;

    size_t i;
    for (i = 0; i < transitionCount(); ++i) {
      st = getTrans(i);
      tt = st.getType();
      w = tt->applyProb(curState,pars);
      time += st.atTime();
      curState += st.getTrans();
      out << prefix << " "
          << name << " " 
          << id << " " << parent << " " << w << " "
          << setw(12) << time << " " 
          << curState << endl;
    }
  }

}


