#ifndef __BRANCHSTATE_H__
#define __BRANCHSTATE_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "MCTraj.h"

namespace MCTraj {
  class BranchState;

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

      int id;
      int old_color;
      int new_color;
  };

  class BranchStates {
    public:
      BranchStates() {}

      BranchStates(size_t n, int col = 0) 
        : colors(n,col), isAlive(n,0)
      {
        alive.reserve(n);
      }

      BranchStates(const BranchStates& bs) 
        : colors(bs.colors), isAlive(bs.isAlive), 
          alive(bs.alive)
      {}

      virtual ~BranchStates() {}

      inline void clear() {
        colors.clear();
        isAlive.clear();
        alive.clear();
      }

      inline void resize(size_t n) {
        clear();
        colors.resize(n,-1);
        isAlive.resize(n,false);
        alive.reserve(n);
      }

      inline BranchStates& operator=(const BranchStates& bs) { 
        colors = bs.colors;
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

      inline int random_alive(gsl_rng* rng) const {
        int id = gsl_rng_uniform_int(rng,alive.size());
        return alive[id];
      }
      int random_color(gsl_rng* rng, int col) const;
      inline size_t num_alive() const { return alive.size(); }
      inline int getCol(size_t i) const { return colors[i]; }

      string to_json() const;

      vector<int> colors;
      vector<bool> isAlive;
      vector<int> alive;
  };
}

#endif
