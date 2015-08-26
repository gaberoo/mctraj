#ifndef __BRANCHSTATE_H__
#define __BRANCHSTATE_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
using namespace std;

#include <rng/Rng.h>

#include "MCTraj.h"

namespace MCTraj {
  class BranchState;

  /**************************************************************************/

  class BranchStateChange {
    friend class BranchState;

    public:
      BranchStateChange() {}

      BranchStateChange(int id, int oc, int nc) 
        : id(id), old_color(oc), new_color(nc) 
      {}

      BranchStateChange(int id, vector<int> c)
        : id(id), change(c)
      {}

      BranchStateChange(const BranchStateChange& bsc) 
        : id(bsc.id), old_color(bsc.old_color), new_color(bsc.new_color),
          change(bsc.change)
      {}

      virtual ~BranchStateChange() {}

      inline BranchStateChange& operator=(const BranchStateChange& bsc) {
        id = bsc.id;
        old_color = bsc.old_color;
        new_color = bsc.new_color;
        change = bsc.change;
        return *this;
      }

      inline void operator=(const vector<int>& v) { change = v; }
      inline void setChange(const int* x, size_t n) { change.assign(x,x+n); }

      template<typename T>
      void json(rapidjson::Writer<T>& json_w) const {
        json_w.StartObject(); {
          json_w.String("id"); json_w.Int(id);
          json_w.String("old"); json_w.Int(old_color);
          json_w.String("new"); json_w.Int(new_color);
          json_w.String("change"); json_w.StartArray(); {
            for (size_t i = 0; i < change.size(); ++i) {
              json_w.Int(change[i]);
            }
          } json_w.EndArray();
        } json_w.EndObject();
      }

      int id;
      int old_color;
      int new_color;
      vector<int> change;
  };

  /**************************************************************************/

  class BranchState : public vector<int> {
    public:
      BranchState() {}
      BranchState(size_t n) : vector<int>(n,0) {}
      virtual ~BranchState() {}

      inline int cnt() const { return std::accumulate(begin(),end(),0); }
  };

  class BranchStates {
    public:
      BranchStates() {}

      BranchStates(size_t n, int col = 0, int nStates = 1) 
        : colors(n,col), states(n,nStates), isAlive(n,0)
      {
        alive.reserve(n);
      }

      BranchStates(const BranchStates& bs) 
        : colors(bs.colors), states(bs.states), 
          isAlive(bs.isAlive), alive(bs.alive)
      {}

      virtual ~BranchStates() {}

      inline void clear() {
        colors.clear();
        states.clear();
        isAlive.clear();
        alive.clear();
      }

      inline void resize(size_t n, size_t nStates = 1) {
        clear();
        colors.resize(n,-1);
        isAlive.resize(n,false);
        alive.reserve(n);
        states.resize(n);
        for (size_t i = 0; i < n; ++i) states[i].resize(nStates);
      }

      inline BranchStates& operator=(const BranchStates& bs) { 
        colors = bs.colors;
        states = bs.states;
        isAlive = bs.isAlive;
        alive = bs.alive;
        return *this;
      }

      BranchStates& operator+=(const BranchStateChange& bsc);
      BranchStates& operator+=(const vector<BranchStateChange>& bsc_vec);

      inline BranchStates& operator-=(const BranchStateChange& bsc) {
        if (bsc.id < (int) colors.size()) colors[bsc.id] = bsc.old_color;
        return *this;
      }

      inline void wake(int i) {
        if (! isAlive[i]) {
          alive.push_back(i);
          isAlive[i] = true;
        }
      }

      inline void kill(int i) {
        if (isAlive[i]) {
          vector<int>::iterator it = find(alive.begin(),alive.end(),i);
          alive.erase(it);
          isAlive[i] = false;
        }
      }

      int countCol(int col) const;
      inline int getCol(size_t i) const { return colors[i]; }
      inline void setCol(size_t i, int c) { colors[i] = c; }
      void aliveCol(int c, vector<int>& alive_col) const;

      inline int getAlive(size_t i) const { return alive[i]; }
      inline size_t nAlive() const { return alive.size(); }
      inline bool awake(size_t id) const { return isAlive[id]; }

      inline int random_alive(rng::RngStream* rng) const {
        int id; 
        rng->uniform_int(1,&id,0,alive.size());
        return alive[id];
      }

      int random_color(rng::RngStream* rng, int col) const;
      inline size_t num_alive() const { return alive.size(); }

      void colWeight(vector<double>& w, int col, bool alive = false) const;
      inline void add(size_t i, size_t c) { states[i][c]++; }
      inline void rem(size_t i, size_t c) { states[i][c]--; }
      inline int state(size_t i, size_t c) const { return states[i][c]; }
      inline int all(size_t i) const { return states[i].cnt(); }

      /* JSON ***************************************************************/

      template<typename T> 
      void json(rapidjson::Writer<T>& json_w, bool aliveOnly = true) const {
        json_w.StartArray();
        if (aliveOnly) {
          for (size_t i = 0; i < alive.size(); ++i) {
            json_w.StartObject(); {
              json_w.String("id");    json_w.Int(alive[i]); 
              json_w.String("color"); json_w.Int(colors[alive[i]]);
              json_w.String("states"); json_w.StartArray(); {
                for (size_t k = 0; k < states[i].size(); ++k) {
                  json_w.Int(states[alive[i]][k]);
                }
              } json_w.EndArray();
            } json_w.EndObject();
          }
        } else {
          for (size_t i = 0; i < colors.size(); ++i) {
            json_w.StartObject(); {
              json_w.String("id");    json_w.Int(i); 
              json_w.String("alive"); json_w.Int(isAlive[i]); 
              json_w.String("color"); json_w.Int(colors[i]);
              json_w.String("states"); json_w.StartArray(); {
                for (size_t k = 0; k < states[i].size(); ++k) {
                  json_w.Int(states[i][k]);
                }
              } json_w.EndArray();
            } json_w.EndObject();
          }
        }
        json_w.EndArray();
      }

      string to_json(bool aliveOnly = true) const;

    public:
      vector<int> colors;
      vector< BranchState > states;

      vector<bool> isAlive;
      vector<int> alive;

      // TODO: keep track of hidden child branches
      // vector< vector<int> > children;
  };
}

#endif
