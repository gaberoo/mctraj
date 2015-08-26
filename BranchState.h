#ifndef __BRANCHSTATE_H__
#define __BRANCHSTATE_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
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

      BranchStateChange(const BranchStateChange& bsc) 
        : id(bsc.id), old_color(bsc.old_color), new_color(bsc.new_color)
      {}

      virtual ~BranchStateChange() {}

      inline BranchStateChange& operator=(const BranchStateChange& bsc) {
        id = bsc.id;
        old_color = bsc.old_color;
        new_color = bsc.new_color;
        return *this;
      }

      template<typename T>
      void json(rapidjson::Writer<T>& json_w) const {
        json_w.StartObject(); {
          json_w.String("id"); json_w.Int(id);
          json_w.String("old"); json_w.Int(old_color);
          json_w.String("new"); json_w.Int(new_color);
        } json_w.EndObject();
      }

      int id;
      int old_color;
      int new_color;
  };

  /**************************************************************************/

  class BranchState : public vector<int> {
    public:
      BranchState() {}
      BranchState(size_t n) : vector<int>(n,0) {}
      virtual ~BranchState() {}
  };

  class BranchStates {
    public:
      BranchStates() {}

      BranchStates(size_t n, int col = 0) 
        : colors(n,col), states(n), isAlive(n,0) 
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

      inline void resize(size_t n) {
        clear();
        colors.resize(n,-1);
        states.resize(n);
        isAlive.resize(n,false);
        alive.reserve(n);
      }

      inline BranchStates& operator=(const BranchStates& bs) { 
        colors = bs.colors;
        states = bs.states;
        isAlive = bs.isAlive;
        alive = bs.alive;
        return *this;
      }

      inline BranchStates& operator+=(const BranchStateChange& bsc) {
        if (bsc.id < (int) colors.size()) colors[bsc.id] = bsc.new_color;
        return *this;
      }
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
      inline int random_alive(rng::RngStream* rng) const {
        int id; 
        rng->uniform_int(1,&id,0,alive.size());
        return alive[id];
      }
      int random_color(rng::RngStream* rng, int col) const;
      inline size_t num_alive() const { return alive.size(); }


      /* JSON ***************************************************************/

      template<typename T> 
      void json(rapidjson::Writer<T>& json_w, bool aliveOnly = true) const {
        json_w.StartArray();
        if (aliveOnly) {
          for (size_t i = 0; i < alive.size(); ++i) {
            json_w.StartObject(); {
              json_w.String("id");    json_w.Int(alive[i]); 
              json_w.String("color"); json_w.Int(colors[alive[i]]);
            } json_w.EndObject();
          }
        } else {
          for (size_t i = 0; i < colors.size(); ++i) {
            json_w.StartObject(); {
              json_w.String("id");    json_w.Int(i); 
              json_w.String("color"); json_w.Int(colors[i]);
            } json_w.EndObject();
          }
        }
        json_w.EndArray();
      }

      string to_json() const;


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
