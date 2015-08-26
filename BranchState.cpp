#include "BranchState.h"

namespace MCTraj {
  BranchStates& BranchStates::operator+=(const BranchStateChange& bsc) {
    if (bsc.id >= 0 && bsc.id < (int) colors.size()) {
      colors[bsc.id] = bsc.new_color;
      if (bsc.new_color < 0) kill(bsc.id);
      else wake(bsc.id);
      if (bsc.change.size() >= states[bsc.id].size()) {
        for (size_t j = 0; j < states[bsc.id].size(); ++j) {
          states[bsc.id][j] += bsc.change[bsc.id];
        }
      }
    }
    return *this;
  }

  BranchStates& BranchStates::operator+=(const vector<BranchStateChange>& bsc) {
    for (size_t i(0); i < bsc.size(); ++i) *this += bsc[i];
    return *this;
  }

  void BranchStates::aliveCol(int col, vector<int>& alive_col) const {
    alive_col.clear();
    for (size_t i = 0; i < alive.size(); ++i) {
      if (colors[alive[i]] == col) alive_col.push_back(alive[i]);
      // cerr << alive[i] << "(" << colors[alive[i]] << "),";
    }
    // cerr << endl;
  }

  void BranchStates::colWeight(vector<double>& w, int col, bool aliveOnly) const {
    w.assign(states.size(),0.0);
    if (aliveOnly) {
      for (size_t i = 0; i < alive.size(); ++i) {
        w[alive[i]] = states[alive[i]].at(col);
      }
      for (size_t i = 1; i < w.size(); ++i) w[i] += w[i-1];
    } else {
      w[0] = states[0].at(col);
      for (size_t i = 1; i < states.size(); ++i) {
        w[i] = w[i-1] + states[i].at(col);
      }
    }
  }

  int BranchStates::random_color(rng::RngStream* rng, int col) const {
    vector<int> col_alive;
    aliveCol(col,col_alive);
    if (col_alive.size() > 0) {
      int r;
      rng->uniform_int(1,&r,0,col_alive.size());
      // cerr << r << " " << col_alive.size() << endl;
      return col_alive[r];
    } else {
      return -1;
    }
  }

  int BranchStates::countCol(int col) const {
    int cnt = 0;
    for (size_t i(0); i < alive.size(); ++i) {
      if (colors[alive[i]] == col) ++cnt;
    }
    return cnt;
  }

  string BranchStates::to_json(bool aliveOnly) const {
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
    json(json_w,aliveOnly);
    return buf.GetString();
  }
}


