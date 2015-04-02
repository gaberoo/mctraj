#include <iostream>
using namespace std;

#include "BranchState.h"
#include <rng/MKLStream.h>

int main() {

  int n = 10;
  MCTraj::BranchStates bs(n);

  bs.wake(4); bs.setCol(4,1);
  bs.wake(6); bs.setCol(6,1);
  bs.wake(8); bs.setCol(8,2);
  bs.wake(9); bs.setCol(9,2);

  rng::RngStream* rng = new rng::MKLStream;
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

