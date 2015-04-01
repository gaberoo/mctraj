#ifndef __TRAJECTORY_H__
#define __TRAJECTORY_H__

#include <iomanip>
#include <vector>
#include <sstream>
using namespace std;

#include <gsl/gsl_rng.h>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

#include "MCTraj.h"
#include "EpiState.h"
#include "TransitionType.h"
#include "StateTransition.h"
#include "Model.h"
#include "TreeNode.h"

namespace MCTraj {
  class Trajectory {
    public:
      Trajectory() 
        : time(0.0), last_event_time(0.0), model(NULL) 
      {}

      Trajectory(const Model* m) 
        : model(m)
      {
        transRates.resize(m->ntrans());
      }

      Trajectory(size_t nstates, const Model* m) 
        : time(0.0), last_event_time(0.0), initialState(nstates), 
          curState(nstates), prob(0.0)
      {
        transRates.resize(m->ntrans());
      }

      Trajectory(const EpiState& es, const Model* m) 
        : time(0.0), last_event_time(0.0), initialState(es), 
          curState(es), model(m), prob(0.0)
      {
        transRates.resize(m->ntrans());
      }

      Trajectory(const Trajectory& T)
        : time(T.time), last_event_time(T.last_event_time),
          initialState(T.initialState), curState(T.curState), 
          transitions(T.transitions),
          model(T.model), transRates(T.transRates),
          prob(T.prob)
      {
        if (model != NULL) transRates.resize(model->ntrans());
      }

      Trajectory(size_t n, istream* filehandle);

      virtual ~Trajectory() {}

      void addTransition(const StateTransition& trans);

      double getTime() const { return time; }

      size_t getState(size_t i) const { return curState[i]; }
      const EpiState& getState() const { return curState; }
      void setState(const EpiState& es) { curState = es; }
      const EpiState& getInitial() const { return initialState; }

      double lastEventTime() const { return last_event_time; }

      inline size_t chooseTransition(gsl_rng* rng) const {
        return model->chooseTransition(rng,transRates);
      }

      int step(double maxTime, const void* pars, gsl_rng* rng, 
               bool noTree = false, bool adjZero = true);
      double force(double nextTime, int nextEvent, const vector<int>& ids, gsl_rng* rng, const void* pars);
      int simulateTrajectory(double endTime, const void* pars, gsl_rng* rng);
      EpiState initState() const { return initialState; }

      size_t transitionCount() const { return transitions.size(); }
      size_t nTrans() const { return transitions.size(); }

      StateTransition& getTrans(size_t i) { return transitions.at(i); }
      const StateTransition& getTrans(size_t i) const { return transitions.at(i); }
      const StateTransition& operator[](size_t i) const { return transitions[i]; }

      const StateTransition& backTrans() const { return transitions.back(); }

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
      inline void resetProb() { prob = 0.0; }

      void printFromLast(size_t last = 0) const;

      ostream& printFromFirst(ostream& out = cout) const;
      ostream& printBranches(ostream& out = cout) const;
      ostream& printTxt(ostream& out = cout) const;

      void toTree(gsl_rng* rng, vector<TreeNode>& tree, const int lineageStates[]) const;

   protected:
      double time;                              /* current process time */
      double last_event_time;                   /* time of last event */

      EpiState initialState;                    /* initial state of trajectory */
      EpiState curState;                        /* current state of the trajectory */

      vector<StateTransition> transitions;      /* transitions */

      const Model* model;
      vector<double> transRates;                /* rates of transition */

      double prob;                              /* probability of trajectory */
  };
}

#endif
