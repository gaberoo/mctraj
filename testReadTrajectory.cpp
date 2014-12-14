#include <iostream>
#include <fstream>
using namespace std;

#include "Trajectory.h"
using namespace MCTraj;


int main(int argc, char** argv) {
  size_t n = 2;
  ifstream fh(argv[1]);
  Trajectory(n,&fh);

  return 0;
}



