#ifndef __TRAJECTORY_H__
#define __TRAJECTORY_H__

#include <iomanip>
#include <vector>
#include <sstream>
#include <set>
#include <numeric>
using namespace std;

#include <rng/Rng.h>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

#include "MCTraj.h"
#include "EpiState.h"
#include "TransitionType.h"
#include "StateTransition.h"
#include "Model.h"
#include "TreeNode.h"

#include "ascii.h"

namespace MCTraj {
  class Trajectory {
    public:
      Trajectory() 
        : time(0.0), last_event_time(0.0), model(NULL),
          store_trans(true)
      {}

      Trajectory(const Model* m) 
        : model(m), store_trans(true)
      {
        transRates.resize(m->nTransRates());
        trueRates.resize(m->nTransRates());
      }

      Trajectory(size_t nstates, const Model* m) 
        : time(0.0), last_event_time(0.0), initialState(nstates), 
          curState(nstates), prob(0.0), store_trans(true)
      {
        transRates.resize(m->nTransRates());
        trueRates.resize(m->nTransRates());
      }

      Trajectory(const EpiState& es, const Model* m) 
        : time(0.0), last_event_time(0.0), initialState(es), 
          curState(es), model(m), prob(0.0), store_trans(true)
      {
        transRates.resize(m->nTransRates());
        trueRates.resize(m->nTransRates());
      }

      Trajectory(const Trajectory& T)
        : time(T.time), last_event_time(T.last_event_time),
          initialState(T.initialState), curState(T.curState), 
          transitions(T.transitions),
          model(T.model), 
          transRates(T.transRates), trueRates(T.trueRates),
          prob(T.prob), store_trans(T.store_trans)
      {
//        if (model != NULL) {
//          transRates.resize(model->nTransRates());
//          trueRates.resize(model->nTransRates());
//        }
      }

      Trajectory(size_t n, istream* filehandle);

      virtual ~Trajectory() {}

      void addTransition(const StateTransition& trans);

      inline double getTime() const { return time; }

      inline size_t getState(size_t i) const { return curState[i]; }
      inline const EpiState& getState() const { return curState; }
      inline void setState(const EpiState& es) { curState = es; }
      inline const EpiState& getInitial() const { return initialState; }

      double lastEventTime() const { return last_event_time; }

      inline size_t chooseTransition(rng::RngStream* rng) const {
        return model->chooseTransition(rng,transRates);
      }

      int step(double maxTime, const void* pars, rng::RngStream* rng, 
               bool adjZero = true);
      double force(double nextTime, int nextEvent, const vector<int>& ids, rng::RngStream* rng, const void* pars);
      int simulateTrajectory(double endTime, const void* pars, rng::RngStream* rng, bool adjZero = true);
      EpiState initState() const { return initialState; }

      inline size_t transitionCount() const { return transitions.size(); }
      inline size_t nTrans() const { return transitions.size(); }

      inline StateTransition& getTrans(size_t i) { return transitions.at(i); }
      inline const StateTransition& getTrans(size_t i) const { return transitions.at(i); }
      inline const StateTransition& backTrans() const { return transitions.back(); }
      inline const StateTransition& operator[](size_t i) const { return transitions[i]; }

      virtual int lastType(int type) const; 

      friend ostream& operator<<(ostream& out, const Trajectory& traj);

      Trajectory& copyState(const Trajectory& T);
      Trajectory& operator=(const Trajectory& T);
      Trajectory& operator+=(const Trajectory& T);

      double allRates(const EpiState& es, void* pars) const;
      double trajProb(const void* pars, int last = -1) const;

      inline double getProb() const { return exp(prob); }
      inline void setProb(double p) { prob = log(p); }
      inline void updateProb(double p) { prob += log(p); }
      inline void updateLogProb(double logp) { prob += logp; }
      inline void resetProb() { prob = 0.0; }

      void printFromLast(size_t last = 0) const;

      ostream& printFromFirst(ostream& out = cout) const;
      ostream& printBranches(ostream& out = cout) const;
      ostream& printTxt(ostream& out = cout) const;

      const Model* getModel() const { return model; }
      inline void storeTrans(bool x) { store_trans = x; }

      void json(rapidjson::Writer<rapidjson::StringBuffer>& w) const;
      string to_json() const;

   protected:
      double time;                         // current process time
      double last_event_time;              // time of last event

      EpiState initialState;               // initial state of trajectory
      EpiState curState;                   // current state of the trajectory
      vector<StateTransition> transitions; // transitions
      const Model* model;                  // dynamic model

      vector<double> transRates;           // simulated rates of transition
      vector<double> trueRates;            // true rates of transition
      double prob;                         // probability of trajectory
      bool store_trans;                    // store simulated trajectory
  };
}

#endif
