#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#ifdef MKLRNG
#include <rng/MKLRng.h>
#else
#include <rng/GSLRng.h>
#endif

#include "pfLik.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"
using namespace MCTraj;

#include <tclap/CmdLine.h>

int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Particle filter approximation for marginal tree likelihood", ' ', "0.9");

  TCLAP::ValueArg<double> N("N","popSize","Total population size",true,100.0,"double",cmd);
  TCLAP::ValueArg<double> beta("b","beta","Transmission rate",true,1.0,"double",cmd);
  TCLAP::ValueArg<double> mu("u","mu","Recovery rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> psi("s","psi","Sequential sampling rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> rho("o","rho","Homochroneous sampling rate",true,0.5,"double",cmd);
  TCLAP::ValueArg<double> gamma("g","gamma","transition rate (for SEIR)",false,0.1,"double",cmd);

  TCLAP::ValueArg<int> numParticles("n","nparticles","Number of particles",false,100,"int",cmd);
  TCLAP::ValueArg<int> reps("r","nreps","Number of repetitions",false,1,"int",cmd);
  TCLAP::ValueArg<int> seed("S","seed","Random number seed",false,-1,"int",cmd);
  TCLAP::ValueArg<int> type("T","model","Model type",false,1,"int",cmd);
  TCLAP::ValueArg<double> alpha("a","alpha","Importance damping",false,10.0,"double",cmd);
  TCLAP::ValueArg<double> stepSize("t","stepSize","Time increments for simulation",false,INFINITY,"double",cmd);

  TCLAP::ValueArg<string> _branchfn("B","branchFn","Filename for branch colors",false,"","string",cmd);
  TCLAP::ValueArg<string> _trajfn("C","trajFn","Filename for trajectory",false,"","string",cmd);

  TCLAP::ValueArg<int> skip("x","skip","Skip lines of times files",false,0,"string",cmd);
  TCLAP::ValueArg<double> filterTime("f","filter","Min time between filters",false,0.0,"double",cmd);

  TCLAP::SwitchArg _printTraj("O","printTraj","Output trajectory",cmd,false);
  TCLAP::SwitchArg printParticles("P","printParticles","Output particles",cmd,false);
  TCLAP::SwitchArg _fullTree("F","fullTree","Full tree",cmd,false);
  TCLAP::SwitchArg adjZero("z","adjZero","Adjust zero-weight trajectory",cmd,false);
  TCLAP::SwitchArg history("H","history","Store trajectories",cmd,false);
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
  pars.N    = N.getValue();
  pars.beta = beta.getValue();
  pars.mu   = mu.getValue();
  pars.psi  = psi.getValue();
  pars.rho  = rho.getValue();

  SEISModel::EpiPars seis_pars;

  PFPars pf_pars;
  pf_pars.num_particles = numParticles.getValue();
  pf_pars.print_particles = printParticles.getValue();
  pf_pars.skip = skip.getValue();
  pf_pars.vflag = vflag.getValue();
  pf_pars.filter_time = filterTime.getValue();
  pf_pars.reps = reps.getValue();
  pf_pars.model_type = type.getValue();
  pf_pars.print_traj = _printTraj.getValue();
  pf_pars.full_tree = _fullTree.getValue();
  pf_pars.seed = (seed.getValue() > 0) ? seed.getValue() : time(NULL);
  pf_pars.adj_zero = adjZero.getValue();
  pf_pars.history = history.getValue();
  pf_pars.step_size = stepSize.getValue();

  string branch_file = _branchfn.getValue();
  string traj_file = _trajfn.getValue();

  // =========================================================================
  
  if (vflag.getValue() > 0) {
    cerr << "# seed        = " << pf_pars.seed << endl;
    cerr << "# adjust zero = " << pf_pars.adj_zero << endl;
  }

  // Setup random number generators

  int max_threads = 1;
#if defined(_OPENMP)
  max_threads = omp_get_max_threads();
#endif

  rng::Rng* rng = NULL;
#ifdef MKLRNG
  rng = new rng::MKLRng;
#else
  rng = new rng::GSLRng;
#endif
  rng->set_seed(pf_pars.seed);
  rng->alloc(max_threads);

  EpiState* es = NULL;
  Model* mpt = NULL;

  switch (pf_pars.model_type) {
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
      seis_pars.gamma = gamma.getValue();
      seis_pars.alpha = alpha.getValue();
      seis_pars.tree = &tree;
      mpt = new SEIS(&seis_pars);
      es = new EpiState(SEISModel::nstates);
      (*es)[0] = ((int) seis_pars.N)-1;
      (*es)[1] = 1;
      (*es)[2] = 0;
      (*es)[3] = 1;
      (*es)[4] = 0;
      if (pf_pars.full_tree) {
        mpt->sim_event(0) = 0;
        mpt->sim_event(1) = 0;
      }
      break;
  }

  double lik = 0.0;

  tree.reverse();
  es->init_branches(tree.max_id()+1,2);
  es->branches.wake(0);
  es->branches.setCol(0,0);
  es->branches.add(0,0);
  // cout << es->to_json() << endl;

  Trajectory* traj = NULL;
  if (pf_pars.print_traj) traj = new Trajectory(*es,mpt);

  ofstream* branch_out = NULL;
  ofstream* traj_out = NULL;

  if (branch_file != "") branch_out = new ofstream(branch_file.c_str());
  if (traj_file != "") traj_out = new ofstream(traj_file.c_str());

  for (int r = 0; r < pf_pars.reps; ++r) {
    lik = pfLik(mpt,*es,tree,pf_pars,rng,traj);
    cout << lik << endl;
    if (branch_out != NULL) traj->printBranches(*branch_out) << endl;
    if (traj_out != NULL) traj->printFromFirst(*traj_out) << endl;
  }

  if (branch_out != NULL) { branch_out->close(); delete branch_out; }
  if (traj_out != NULL) { traj_out->close(); delete traj_out; }

  // if (pf_pars.print_traj) cout << *traj << endl;

  delete traj;
  delete mpt;
  delete es;

  return 0;
}
