#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <gsl/gsl_rng.h>

#include "pfLik.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"

using namespace MCTraj;

int main(int argc, char** argv) {
  SISModel::EpiPars pars = { 100.0, 1.0, 0.1, 0.1, 1.0 };
  SEISModel::EpiPars seis_pars;
  double _gamma = 0.0;

  size_t num_particles = 1;
  int reps = 1;
  int modelType = 1;
  int printTraj = 0;
  int printParticles = 0;
  int fullTree = 0;
  long unsigned seed = time(NULL);
  int skip = 1;

  string branch_file = "";
  string traj_file = "";

  int vflag = 0;
  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"N:b:u:s:o:g:vn:r:ROPS:T:B:C:Fx:")) != -1) {
    switch (c) {
      case 'N': pars.N = atof(optarg); break;
      case 'b': pars.beta = atof(optarg); break;
      case 'u': pars.mu = atof(optarg); break;
      case 's': pars.psi = atof(optarg); break;
      case 'o': pars.rho = atof(optarg); break;
      case 'v': ++vflag; break;
      case 'n': num_particles = atoi(optarg); break;
      case 'r': reps = atoi(optarg); break;
      case 'O': printTraj = ! printTraj; break;
      case 'P': printParticles = ! printParticles; break;
      case 'S': seed = atoi(optarg); break;
      case 'T': modelType = atoi(optarg); break;
      case 'g': _gamma = atof(optarg); break;
      case 'B': branch_file = optarg; break;
      case 'C': traj_file = optarg; break;
      case 'F': fullTree = ! fullTree; break;
      case 'x': skip = atoi(optarg); break;
    }
  }

  if (optind == argc) {
    cout << "Please supply a tree file." << endl;
    return 0;
  }
  Tree tree(*(argv+optind));

  // Setup random number generators
  int max_threads = omp_get_max_threads();
  gsl_rng** rng = (gsl_rng**) malloc(max_threads*sizeof(gsl_rng*));
  for (int i = 0; i < max_threads; ++i) {
    rng[i] = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng[i],seed+i);
  }

  EpiState* es = NULL;
  Model* mpt = NULL;

  switch (modelType) {
    case 1:
    default:
      mpt = new SIS(&pars);
      es = new EpiState(SISModel::nstates);
      (*es)[0] = (int) pars.N;
      break;
    case 2:
      mpt = new SIR(&pars);
      es = new EpiState(SISModel::nstates);
      (*es)[0] = (int) pars.N;
      break;
    case 3:
      seis_pars.N = pars.N;
      seis_pars.beta = pars.beta;
      seis_pars.mu = pars.mu;
      seis_pars.psi = pars.psi;
      seis_pars.rho = pars.rho;
      seis_pars.gamma = _gamma;
      mpt = new SEIS(&seis_pars);
      es = new EpiState(SEISModel::nstates);
      (*es)[0] = ((int) seis_pars.N)-1;
      (*es)[1] = 1;
      (*es)[2] = 0;
      (*es)[3] = 1;
      (*es)[4] = 0;
      if (fullTree) {
        mpt->sim_event(0) = 0;
        mpt->sim_event(1) = 0;
      }
      break;
  }

  double lik = 0.0;

  tree.reverse();
  es->init_branches(tree.max_id()+1);
  es->branches.wake(0);
  es->branches.setCol(0,0);
  // cout << es->to_json() << endl;

  Trajectory* traj = NULL;
  if (printTraj) traj = new Trajectory(*es,mpt);

  ofstream* branch_out = NULL;
  ofstream* traj_out = NULL;

  if (branch_file != "") branch_out = new ofstream(branch_file.c_str());
  if (traj_file != "") traj_out = new ofstream(traj_file.c_str());

  for (int r = 0; r < reps; ++r) {
    lik = pfLik(mpt,*es,tree,num_particles,rng,vflag,traj,skip,printParticles);
    cout << lik << endl;
    if (branch_out != NULL) traj->printBranches(*branch_out) << endl;
    if (traj_out != NULL) traj->printFromFirst(*traj_out) << endl;
  }

  if (branch_out != NULL) { branch_out->close(); delete branch_out; }
  if (traj_out != NULL) { traj_out->close(); delete traj_out; }

  if (printTraj) cout << *traj << endl;

  delete traj;
  for (int i = 0; i < max_threads; ++i) gsl_rng_free(rng[i]);
  free(rng);
  delete mpt;
  delete es;

  return 0;
}
