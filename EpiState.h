#ifndef __EPISTATE_H__
#define __EPISTATE_H__

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
using namespace std;

#include "MCTraj.h"
#include "BranchState.h"

#include <rapidjson/document.h>

namespace MCTraj {
  class EpiState {
    public:
      EpiState() {}
      EpiState(size_t nstates) : state(nstates), time(0.0) {}
      EpiState(size_t nstates, int x) : state(nstates,x), time(0.0) {}
      EpiState(const vector<int>& es) : state(es), time(0.0) {}
      EpiState(const EpiState& es) 
        : state(es.state), 
          time(es.time),
          branches(es.branches), 
          curBranch(es.curBranch)
      {}
      virtual ~EpiState() {}

      inline size_t numStates() const { return state.size(); }
      inline int incease(size_t i) { return ++state[i]; }
      inline int decrease(size_t i) { return --state[i]; }
      inline int at(size_t i) const { return state.at(i); }
      inline int& operator[](size_t i) { return state[i]; }
      inline int operator[](size_t i) const { return state[i]; }

      EpiState& apply(const vector<int>& change, int dir = 1) {
        if (state.size() <= change.size()) {
          for (size_t i(0); i < state.size(); ++i) {
            state[i] += dir*change[i];
          }
        }
        return *this;
      }

      const vector<int>& s() const { return state; }

      EpiState& operator=(const EpiState& es) { 
        time = es.time;
        state = es.state; 
        branches = es.branches;
        curBranch = es.curBranch;
        return *this; 
      }

      EpiState& operator+=(const vector<int>& change) { return apply(change); }
      EpiState& operator-=(const vector<int>& change) { return apply(change,-1); }

      friend ostream& operator<<(ostream& out, const EpiState& es) {
        for (size_t i(0); i < es.state.size(); ++i) {
          out << setw(4) << es.state[i] << " ";
        }
        return out;
      }

      void init_branches(size_t n) { branches.resize(n); }

      /* JSON ***************************************************************/

      template<typename T> 
      void json(rapidjson::Writer<T>& json_w) const {
        json_w.StartObject(); {
          json_w.String("pop"); json_w.StartArray();
          for (size_t i(0); i < state.size(); ++i) json_w.Int(state[i]);
          json_w.EndArray();
          json_w.String("branches");
          branches.json(json_w);
        } json_w.EndObject();
      }
      string to_json() const;

    protected:
      vector<int> state;

    public:
      double time;
      BranchStates branches;
      vector<int> curBranch;
  };
}

#endif

