#include "pfLik.h"
#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

int main(int argc, char** argv) {
  SISModel::EpiPars pars = { 100.0, 1.0, 0.1, 0.1 };
  size_t num_particles = 1;

  long unsigned int seed = 0;
  int vflag = 0;
  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"N:b:u:s:vn:r:t:S:")) != -1) {
    switch (c) {
      case 'N': pars.N = atof(optarg); break;
      case 'b': pars.beta = atof(optarg); break;
      case 'u': pars.mu = atof(optarg); break;
      case 's': pars.psi = atof(optarg); break;
      case 'v': ++vflag; break;
      case 'n': num_particles = atoi(optarg); break;
      case 'S': seed = atoi(optarg); break;
    }
  }

  size_t nstates = 3;
  EpiState init(nstates);
  init[1] = 0;
  init[0] = (int) (pars.N - init[1]);
  init[2] = 0;

  Tree tree(*(argv+optind));
  tree.reverse();

    // cout << "Starting PF likelihood calculation..." << endl;

  if (seed == 0) seed = time(NULL);
  int max_threads = omp_get_max_threads();
  gsl_rng** rng = (gsl_rng**) malloc(max_threads*sizeof(gsl_rng*));
  for (int i = 0; i < max_threads; ++i) {
    rng[i] = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng[i],seed+i);
  }

  SIS m(&pars);
  Trajectory T(init,&m);
  TrajParticleFilter pf(&m);

  for (size_t i(0); i < num_particles; ++i) {
    char name[32];
    sprintf(name,"P%d",(int) i);
    pf.push_back(TrajParticle(string(name),1.0,T));
    pf[i].setId(i);
  }
  pf.setTree(&tree,1);

  int ret = 1;
  int filterRet = 0;

  while (ret >= 0) {
    ret = pf.stepTree(&pars,rng);

    if (ret >= 0) {
      pf.calcWeights(&pars);
      pf.addTreeEvent(&pars);
      filterRet = pf.filter(rng);

      if (filterRet < 0) {
        cerr << "\033[1;31m";
        cerr << pf[0].getTime() << "/" << pf.maxTime() << ": Particle collapse!";
        cerr << "\033[0m" << endl;
        return -INFINITY;
      }
    }
  }

  cerr << "\033[1;32m";
  cerr << "(" << pars.N << "," << pars.beta << "," << pars.mu << "," << pars.psi << ") ";
  cerr << "LL = " << pf.est() << " / " << pf.est_p();
  cerr << "\033[0m" << endl;

  StateTransition st;
  vector<const TrajParticle*> traj(pf.singleTraj(0));
  // vector<const TrajParticle*>::const_reverse_iterator t;
  const TrajParticle* tp;
  size_t j;
  for (j = traj.size()-1; j > 0; --j) {
    tp = traj.at(j);
    double time = tp->getInit();
    EpiState curState = tp->getInitial();
    cout << "+ ";
    cout << setw(12) << time << " " << curState;
    // if (j > 0) cout << " " << tree.time(j-1) << " " << tree.ttype(j-1);
    int k = traj.size()-j-2;
    if (k >= 0 && k < (int) tree.size()) 
      cout << " " << setw(8) << tree.time(k) << " " << tree.ttype(k);
    else 
      cout << " " << -1 << " " << -1;
    cout << endl;
    for (size_t i = 0; i < tp->transitionCount()-1; ++i) {
      st = tp->getTrans(i);
      time += st.atTime();
      curState += st.getTrans();
      cout << "- ";
      cout << setw(12) << time << " " << curState;
      cout << " " << 0 << " " << st.etype() << endl;
    }
  }

  for (int i = 0; i < max_threads; ++i) gsl_rng_free(rng[i]);
  free(rng);

  return 0;
}
