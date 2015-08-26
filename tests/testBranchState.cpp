#include <iostream>
using namespace std;

#include "BranchState.h"
#include <rng/GSLStream.h>

int main() {

  int n = 10;
  MCTraj::BranchStates bs(n,0,2);

  bs.wake(4); bs.setCol(4,1);
  bs.wake(6); bs.setCol(6,1);
  bs.wake(8); bs.setCol(8,2);
  bs.wake(9); bs.setCol(9,2);

  int i;

  i = 0; while (i++ < 3) bs.add(4,0);
  i = 0; while (i++ < 2) bs.add(4,1);

  bs.add(6,0); bs.add(6,1);
  bs.add(5,0); bs.add(5,1);

  for (int i = 0; i < n; ++i) {
    cout << i << " " << bs.colors[i] << " " << bs.isAlive[i] << endl;
  }
  cout << endl;

  cout << bs.to_json(false) << endl;

  vector<double> w;

  bs.colWeight(w,1,true);
  for (size_t i = 0; i < w.size(); ++i) {
    cout << i << " " << w[i] << endl;
  }

  bs.colWeight(w,1,false);
  for (size_t i = 0; i < w.size(); ++i) {
    cout << i << " " << w[i] << endl;
  }

  rng::RngStream* rng = new rng::GSLStream;
  rng->alloc(time(NULL));

  int r[10];
  rng->uniform_int(10,r,0,2);
  for (int i = 0; i < 10; ++i) cout << r[i] << " ";
  cout << endl;

  int c = bs.random_color(rng,1);
  cout << c << endl;

  int cnt = bs.countCol(1);
  cout << cnt << endl;

  cout << bs.to_json() << endl;

  return 0;
}

