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

  class BranchStates {
    public:
      BranchStates() 
        : size(0), nstates(1)
      {}

      BranchStates(size_t n, int col = 0, int nStates = 1) 
        : size(n), nstates(nStates),
          colors(n,col), v_states(n*(nStates+1),0), 
          isAlive(n,false)
      {
        alive.reserve(n);
      }

      BranchStates(const BranchStates& bs) 
        : size(bs.size), nstates(bs.nstates),
          colors(bs.colors), v_states(bs.v_states), 
          isAlive(bs.isAlive), alive(bs.alive)
      {}

      virtual ~BranchStates() {}

      inline void clear() {
        colors.clear();
        v_states.clear();
        isAlive.clear();
        alive.clear();
      }

      inline void resize(size_t n, size_t nStates = 1) {
        clear();
        colors.resize(n,-1);
        isAlive.resize(n,false);
        alive.reserve(n);
        size = n;
        nstates = nStates;
        v_states.assign(n*(nStates+1),0);
      }

      inline BranchStates& operator=(const BranchStates& bs) { 
        size = bs.size;
        nstates = bs.nstates;
        colors = bs.colors;
        v_states = bs.v_states;
        isAlive = bs.isAlive;
        alive = bs.alive;
        return *this;
      }

      BranchStates& operator+=(const BranchStateChange& bsc);
      BranchStates& operator+=(const vector<BranchStateChange>& bsc_vec);
      // BranchStates& operator-=(const BranchStateChange& bsc);

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

      int colProb(vector<double>& w, int col, bool alive = false) const;
      int colWeight(vector<double>& w, int col, bool alive = false) const;

//      inline const int* states(size_t i) const { return v_states.data() + (i*size); }
//      inline int* states(size_t i) { return v_states.data() + (i*size); }

      inline int state(size_t i, size_t c) const { return v_states.at(i*(nstates+1)+c); }
      inline int& state(size_t i, size_t c) { return v_states.at(i*(nstates+1)+c); }

      inline void add(size_t i, size_t c) { state(i,c)++; }
      inline void rem(size_t i, size_t c) { state(i,c)--; }

      inline int all(size_t i) const {
        vector<int>::const_iterator it = v_states.begin() + i*(nstates+1) + 1;
        return std::accumulate(it,it+nstates,0);
      }

      /* JSON ***************************************************************/

      template<typename T> 
      void json(rapidjson::Writer<T>& json_w, bool aliveOnly = true) const {
        json_w.StartArray();
        if (aliveOnly) {
          for (size_t i = 0; i < alive.size(); ++i) {
            state_json(i,json_w);
          }
        } else {
          for (size_t i = 0; i < colors.size(); ++i) {
            state_json(i,json_w);
          }
        }
        json_w.EndArray();
      }

      template<typename T> 
      void state_json(size_t i, rapidjson::Writer<T>& json_w) const {
        if (i < size) {
          json_w.StartObject(); {
            json_w.String("id");    json_w.Int(i); 
            json_w.String("alive"); json_w.Int(isAlive[i]); 
            json_w.String("color"); json_w.Int(colors[i]);
            json_w.String("states");
            json_w.StartArray();
            for (int k = 0; k <= nstates; ++k) json_w.Int(state(i,k));
            json_w.EndArray();
          } json_w.EndObject();
        } else {
          json_w.String("Err");
        }
      }

      string to_json(bool aliveOnly = true) const;
      string state_to_json(int id) const;

    public:
      int size;
      int nstates;

      vector<int> colors;
      vector<int> v_states;

      vector<bool> isAlive;
      vector<int> alive;

      // vector< BranchState > states;
  };
}

#endif
