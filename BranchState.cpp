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

  int BranchStates::random_color(rng::RngStream* rng, int col) const {
    vector<int> col_alive;
    for (size_t i(0); i < alive.size(); ++i) {
      if (colors[alive[i]] == col) col_alive.push_back(alive[i]);
      // cerr << alive[i] << "(" << colors[alive[i]] << "),";
    }
    // cerr << endl;
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

  string BranchStates::to_json() const {
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
    json(json_w,true);
    return buf.GetString();
  }
}


