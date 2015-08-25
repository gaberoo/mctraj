#include <iostream>
using namespace std;

#include <rng/GSLRng.h>
#include <rng/MKLRng.h>
using namespace rng;

int main() {
  // Rng* rng = new MKLRng;
  Rng* rng = new GSLRng;

  unsigned long seed = time(NULL);
  rng->set_seed(seed);
  rng->alloc(1);

  double r = 0.0;
  rng->uniform(0,1,&r);
  rng->uniform(0,1,&r);
  // cout << seed << " " << r << endl;

//  double w[] = { 0.1, 0.2, 0.2, 0.5 };
//  for (int i = 0; i < 100000; ++i) {
//    cout << (*rng)[0]->discrete(4,w) << endl;
//  }

  double p[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  for (int i = 0; i < 100000; ++i) {
    cout << (*rng)[0]->discrete_x(5,p) << " "
         << (*rng)[0]->pick(p,5) << endl;
  }

  rng->free();
  return 0;
}

