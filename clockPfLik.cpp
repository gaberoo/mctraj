#include "pfLik.h"
#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

int main(int argc, char** argv) {
  SISModel::EpiPars pars = { 100.0, 1.0, 0.1, 0.1 };
  size_t num_particles = 1;

  int vflag = 0;
  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"N:b:u:s:vn:r:t:")) != -1) {
    switch (c) {
      case 'N': pars.N = atof(optarg); break;
      case 'b': pars.beta = atof(optarg); break;
      case 'u': pars.mu = atof(optarg); break;
      case 's': pars.psi = atof(optarg); break;
      case 'v': ++vflag; break;
      case 'n': num_particles = atoi(optarg); break;
    }
  }

  size_t nstates = 2;
  EpiState init(nstates);
  init[1] = 1;
  init[0] = (int) (pars.N - init[1]);

  Tree tree(*(argv+optind));
  tree.reverse();

    // cout << "Starting PF likelihood calculation..." << endl;

  SIS m(&pars);
  Trajectory T(init,&m);
  TrajParticleFilter pf(&m);

  for (size_t i(0); i < num_particles; ++i) {
    char name[32];
    sprintf(name,"P%d",(int) i);
    // cout << "  " << name << endl;
    pf.push_back(TrajParticle(string(name),1.0,T));
    pf[i].setId(i);
  }
  pf.setTree(&tree);

  int ret = 1;
  int filterRet;

  double t;
  double ct[] = { 0.0, 0.0, 0.0, 0.0 };

  unsigned long seed = time(NULL);
  int max_threads = omp_get_max_threads();
  gsl_rng** rng = (gsl_rng**) malloc(max_threads*sizeof(gsl_rng*));
  for (int i = 0; i < max_threads; ++i) {
    rng[i] = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng[i],seed+i);
  }

  while (ret > 0) {
    t = omp_get_wtime();
    ret = pf.stepTree(&pars,rng);
    ct[0] += omp_get_wtime() - t;

    if (ret > 0) {
      t = omp_get_wtime();
      pf.calcWeights(&pars);
      ct[1] += omp_get_wtime() - t;

      t = omp_get_wtime();
      pf.addTreeEvent(&pars);
      ct[2] += omp_get_wtime() - t;

      t = omp_get_wtime();
      filterRet = pf.filter(rng,'p');
      ct[3] += omp_get_wtime() - t;

      if (filterRet < 0) {
        cerr << "\033[1;31m";
        cerr << pf[0].getTime() << "/" << pf.maxTime() << ": Particle collapse!";
        cerr << "\033[0m" << endl;
        return -INFINITY;
      }
    }
  }

  cerr << "Calculation times:" << endl;
  cerr << "  step tree      : " << ct[0] << endl;
  cerr << "  calc weights   : " << ct[1] << endl;
  cerr << "  add tree event : " << ct[2] << endl;
  cerr << "  filter         : " << ct[3] << endl;

  cerr << "\033[1;32m";
  cerr << "(" << pars.N << "," << pars.beta << "," << pars.mu << "," << pars.psi << ") ";
  cerr << "LL = " << pf.est() << " / " << pf.est_p();
  cerr << "\033[0m" << endl;

  for (int i = 0; i < max_threads; ++i) gsl_rng_free(rng[i]);
  free(rng);

  return 0;
}
