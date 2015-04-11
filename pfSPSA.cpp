#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <rng/GSLRng.h>
#include <spsa/spsa.h>

#include "pfLik.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"
using namespace MCTraj;

#include <tclap/CmdLine.h>

typedef struct {
  Model* mpt;
  const EpiState* es;
  const Tree* tree;
  const PFPars* pars;
  rng::Rng* rng;
  const int* lock;
} pf_pars_t;

double pf_lik(const double* state, void* pars) 
{
  pf_pars_t& p = *(pf_pars_t*) pars;

  const Pars* oldPt = p.mpt->getPars();
  const SEISModel::EpiPars* oldPars;

  try {
    oldPars = dynamic_cast<const SEISModel::EpiPars*>(oldPt);
  } catch (exception& e) {
    cerr << "Error casting pointer!" << endl;
    abort();
  }

  SEISModel::EpiPars epi(*oldPars);
  size_t i = 0;
  if (! p.lock[0]) epi.N = exp(state[i++]);
  if (! p.lock[1]) epi.beta = exp(state[i++]);
  if (! p.lock[2]) epi.mu = exp(state[i++]);
  if (! p.lock[3]) epi.psi = exp(state[i++]);
  if (! p.lock[4]) epi.gamma = exp(state[i++]);

  EpiState init(*p.es);
  init[0] = (int) ceil(epi.N) - 1;
  if (init[0] < 0) return -INFINITY;

  p.mpt->setPars(&epi);
  double lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,NULL);
  p.mpt->setPars(oldPt);

  return -lik;
}

int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Particle filter approximation for marginal tree likelihood", ' ', "0.9");

  TCLAP::ValueArg<double> N("N","popSize","Total population size",true,100.0,"double",cmd);
  TCLAP::ValueArg<double> beta("b","beta","Transmission rate",true,1.0,"double",cmd);
  TCLAP::ValueArg<double> mu("u","mu","Recovery rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> psi("s","psi","Sequential sampling rate",false,0.1,"double",cmd);
  TCLAP::ValueArg<double> rho("o","rho","Homochroneous sampling rate",false,0.0,"double",cmd);
  TCLAP::ValueArg<double> gamma("g","gamma","transition rate (for SEIR)",false,0.1,"double",cmd);

  TCLAP::ValueArg<int> numParticles("n","nparticles","Number of particles",false,100,"int",cmd);
  TCLAP::ValueArg<int> reps("r","nreps","Number of repetitions",false,1,"int",cmd);
  TCLAP::ValueArg<int> seed("S","seed","Random number seed",false,-1,"int",cmd);
  TCLAP::ValueArg<int> type("T","model","Model type",false,1,"int",cmd);
  TCLAP::ValueArg<int> skip("x","skip","Skip lines of times files",false,0,"string",cmd);

  TCLAP::ValueArg<double> filterTime("f","filter","Min time between filters",false,0.0,"double",cmd);

  TCLAP::ValueArg<int> steps("K","steps","Number of SPSA iterations",false,100,"int",cmd);
  TCLAP::ValueArg<double> spsa_a("a","ak","Step size for gradient method",false,0.01,"double",cmd);
  TCLAP::ValueArg<double> spsa_c("c","ck","Step size for gradient method",false,1.0,"double",cmd);

  TCLAP::SwitchArg _printTraj("O","printTraj","Output trajectory",cmd,false);
  TCLAP::SwitchArg printParticles("P","printParticles","Output particles",cmd,false);
  TCLAP::SwitchArg _fullTree("F","fullTree","Full tree",cmd,false);
  TCLAP::SwitchArg adjZero("z","adjZero","Adjust zero-weight trajectory",cmd,false);
  TCLAP::SwitchArg calcLik("L","calcLik","Calculate likelihood of each estimate",cmd,false);
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

  // =========================================================================
  
  // if (vflag.getValue() > 0) cerr << "# seed = " << seed << endl;

  // Setup random number generators
  int max_threads = omp_get_max_threads();
  rng::Rng* rng = new rng::GSLRng;
  rng->set_seed(pf_pars.seed);
  rng->alloc(max_threads);

  SEISModel::EpiPars seis_pars;
  seis_pars.N    = fabs(N.getValue());
  seis_pars.beta = fabs(beta.getValue());
  seis_pars.mu   = fabs(mu.getValue());
  seis_pars.psi  = fabs(psi.getValue());
  seis_pars.rho  = fabs(rho.getValue());
  seis_pars.gamma = fabs(gamma.getValue());

  Model* mpt = new SEIS(&seis_pars);
  EpiState* es = new EpiState(SEISModel::nstates);
  (*es)[0] = ((int) seis_pars.N)-1;
  (*es)[1] = 1;
  (*es)[2] = 0;
  (*es)[3] = 1;
  (*es)[4] = 0;

  tree.reverse();
  es->init_branches(tree.max_id()+1);
  es->branches.wake(0);
  es->branches.setCol(0,0);
  // cout << es->to_json() << endl;

  // Start of SPSA ===========================================================

  // PF parameters
  pf_pars_t pfp;
  pfp.mpt = mpt;
  pfp.es = es;
  pfp.tree = &tree;
  pfp.pars = &pf_pars;
  pfp.rng = rng;

  // SPSA parameters
  spsa::pars_t p;
  p.a = spsa_a.getValue();
  p.c = spsa_c.getValue();
  p.alpha = 0.602;
  p.gamma = 0.101;
  p.A = 5.0;
  p.rng = (*rng)[0];
  p.fun = &pf_lik;
  p.ak = 0.0;
  p.ck = 0.0;
  p.pars = &pfp;

  int lock[] = { 0, 0, 0, 0, 0 };

  size_t nvar = 0;
  double theta[5];

  if (N.getValue() <= 0.0)     { lock[0] = 1; } else { theta[nvar++] = log(N.getValue()); }
  if (beta.getValue() <= 0.0)  { lock[1] = 1; } else { theta[nvar++] = log(beta.getValue()); }
  if (mu.getValue() <= 0.0)    { lock[2] = 1; } else { theta[nvar++] = log(mu.getValue()); }
  if (psi.getValue() <= 0.0)   { lock[3] = 1; } else { theta[nvar++] = log(psi.getValue()); }
  if (gamma.getValue() <= 0.0) { lock[4] = 1; } else { theta[nvar++] = log(gamma.getValue()); }

  if (nvar <= 0) {
    cerr << "All variables locked!" << endl;
    return 0;
  }

  pfp.lock = lock;

  double lik = GSL_NAN;

  if (calcLik.getValue()) lik = pf_lik(theta,&pfp);

  spsa::init_pars(&p,nvar);
  // cerr << " Lik = " << pf_lik(theta,&pfp) << endl;

  int ret = 0;
  for (int k = 1; k <= steps.getValue(); ++k) {
    p.ak = p.a/pow(k+1.0+p.A,p.alpha);
    p.ck = p.c/pow(k+1.0,p.gamma);
    ret = spsa::approx_gradient(&p,theta);
    if (ret > 0) {
      for (size_t i = 0; i < p.n; ++i) {
        theta[i] -= p.ak*p.grad[i];
      }
      if (calcLik.getValue()) {
        switch (ret) {
          case 1:
            lik = p.y[1];
            break;
          case 2:
            lik = p.y[0];
            break;
          case 3:
          default:
            lik = pf_lik(theta,&pfp);
            break;
        }
      }
    } else {
      cerr << "ROOOAAAAARRRR -- something went wrong in evaluating the likelihood." << endl;
      return 0;
    }
    cout << k << " " << ret << " " << lik << " ";
    for (size_t i = 0; i < p.n; ++i) cout << exp(theta[i]) << " ";
    cout << p.y[0] << " ";
    for (size_t i = 0; i < p.n; ++i) cout << exp(p.x1[i]) << " ";
    cout << p.y[1] << " ";
    for (size_t i = 0; i < p.n; ++i) cout << exp(p.x2[i]) << " ";
    cout << endl;
  }

  spsa::free_pars(&p);

  // End of SPSA =============================================================

  delete mpt;
  delete es;

  return 0;
}
