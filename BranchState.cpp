#include "BranchState.h"

namespace MCTraj {
  BranchStates& BranchStates::operator+=(const vector<BranchStateChange>& bsc) {
    for (size_t i(0); i < bsc.size(); ++i) {
      if (bsc[i].id >= 0 && bsc[i].id < (int) colors.size()) {
        colors[bsc[i].id] = bsc[i].new_color;
        if (bsc[i].new_color < 0) kill(bsc[i].id);
        else wake(bsc[i].id);
        // cerr << "   " << bsc[i].id << " -> " << colors[bsc[i].id] << endl;
      }
    }
    return *this;
  }

  int BranchStates::random_color(gsl_rng* rng, int col) const {
    vector<int> col_alive;
    for (size_t i(0); i < alive.size(); ++i) {
      if (colors[alive[i]] == col) col_alive.push_back(alive[i]);
      // cerr << alive[i] << "(" << colors[alive[i]] << "),";
    }
    // cerr << endl;
    if (col_alive.size() > 0) {
      int r = gsl_rng_uniform_int(rng,col_alive.size());
      return col_alive[r];
    } else {
      return -1;
    }
  }

  string BranchStates::to_json() const {
    ostringstream out;
    size_t i;
    out << "{";
    for (i = 0; i < alive.size(); ++i) {
      out << alive[i] << ":" << colors[alive[i]];
      if (i < alive.size()-1) out << ",";
    }
    out << "}";
    return out.str();
  }
}


