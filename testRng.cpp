#include <iostream>
using namespace std;

#include <rng/MKLRng.h>
using namespace rng;

int main() {
  Rng* rng = new MKLRng;

  unsigned long seed = time(NULL);
  rng->set_seed(seed);
  rng->alloc(1);

  double r = 0.0;
  rng->uniform(0,1,&r);
  rng->uniform(0,1,&r);
  // cout << seed << " " << r << endl;

  double w[] = { 0.1, 0.2, 0.2, 0.5 };
  for (int i = 0; i < 100000; ++i) {
    cout << (*rng)[0]->discrete(4,w) << endl;
  }

  rng->free();
  return 0;
}

