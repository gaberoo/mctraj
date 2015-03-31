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

#include <tclap/CmdLine.h>

int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Particle filter approximation for marginal tree likelihood", ' ', "0.9");

  TCLAP::ValueArg<double> _N("N","popSize","Total population size",true,100.0,"double",cmd);
  TCLAP::ValueArg<double> _beta("b","beta","Transmission rate",true,1.0,"double",cmd);
  TCLAP::ValueArg<double> _mu("u","mu","Recovery rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> _psi("s","psi","Sequential sampling rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> _rho("o","rho","Homochroneous sampling rate",true,0.5,"double",cmd);
  TCLAP::ValueArg<double> _gamma("g","gamma","transition rate (for SEIR)",false,0.1,"double",cmd);
  TCLAP::ValueArg<int> _numParticles("n","nparticles","Number of particles",false,100,"int",cmd);
  TCLAP::ValueArg<int> _reps("r","nreps","Number of repetitions",false,1,"int",cmd);
  TCLAP::ValueArg<int> _seed("S","seed","Random number seed",false,-1,"int",cmd);
  TCLAP::ValueArg<int> _type("T","model","Model type",false,1,"int",cmd);
  TCLAP::ValueArg<string> _branchfn("B","branchFn","Filename for branch colors",false,"","string",cmd);
  TCLAP::ValueArg<string> _trajfn("C","trajFn","Filename for trajectory",false,"","string",cmd);
  TCLAP::ValueArg<int> _skip("x","skip","Skip lines of times files",false,0,"string",cmd);
  TCLAP::ValueArg<double> filterTime("f","filter","Min time between filters",false,0.0,"double",cmd);

  TCLAP::SwitchArg _printTraj("O","printTraj","Output trajectory",cmd,false);
  TCLAP::SwitchArg _printParticles("P","printParticles","Output particles",cmd,false);
  TCLAP::SwitchArg _fullTree("F","fullTree","Full tree",cmd,false);
  TCLAP::MultiSwitchArg vflag("v","verbose","Increase verbosity",cmd);

  TCLAP::UnlabeledMultiArg<string> multi("files","Trees",true,"Input Tree files",false);
  cmd.add(multi);

  try {
    cmd.parse(argc,argv);
  } catch (TCLAP::ArgException &e) { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

  if (optind == argc) {
    cout << "Please supply a tree file." << endl;
    return 0;
  }

  vector<string> fileNames = multi.getValue();
  Tree tree(fileNames.front().c_str());

  // =========================================================================

  SISModel::EpiPars pars;
  pars.N    = _N.getValue();
  pars.beta = _beta.getValue();
  pars.mu   = _mu.getValue();
  pars.psi  = _psi.getValue();
  pars.rho  = _rho.getValue();

  SEISModel::EpiPars seis_pars;

  size_t num_particles = _numParticles.getValue();
  int reps = _reps.getValue();
  int modelType = _type.getValue();
  int printTraj = _printTraj.getValue();
  int printParticles = _printParticles.getValue();
  int fullTree = _fullTree.getValue();
  long unsigned seed = (_seed.getValue() > 0) ? _seed.getValue() : time(NULL);
  int skip = _skip.getValue();

  string branch_file = _branchfn.getValue();
  string traj_file = _trajfn.getValue();

  // =========================================================================
  
  // if (vflag.getValue() > 0) cerr << "# seed = " << seed << endl;

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
      seis_pars.gamma = _gamma.getValue();
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
    lik = pfLik(mpt,*es,tree,num_particles,rng,filterTime.getValue(),
                vflag.getValue(),traj,skip,printParticles);
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
