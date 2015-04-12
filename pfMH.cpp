#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <tclap/CmdLine.h>
#include <rapidjson/document.h>

#include <rng/GSLRng.h>
#include <spsa/spsa.h>
#include <cdream/dream.h>

#include "pfLik.h"
#include "models/SIS.h"
#include "models/SIR.h"
#include "models/SEIS.h"
using namespace MCTraj;

typedef struct {
  Model* mpt;
  const Tree* tree;
  const PFPars* pars;
  const EpiState* es;
  rng::Rng* rng;
  const dream_pars* mpar;
} pf_pars_t;

double pf_lik(const double* state, const void* pars) 
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
  epi.N     = (p.mpar->scale[0] == 'l') ? exp(state[0]) : state[0];
  epi.beta  = (p.mpar->scale[1] == 'l') ? exp(state[1]) : state[1];
  epi.mu    = (p.mpar->scale[2] == 'l') ? exp(state[2]) : state[2];
  epi.psi   = (p.mpar->scale[3] == 'l') ? exp(state[3]) : state[3];
  epi.gamma = (p.mpar->scale[4] == 'l') ? exp(state[4]) : state[4];

  EpiState init(*p.es);
  init[0] = (int) ceil(epi.N) - 1;
  if (init[0] < 0) return -INFINITY;

  if (p.pars->vflag) cerr << epi.to_json() << endl;
  p.mpt->setPars(&epi);

  double lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,NULL);

  p.mpt->setPars(oldPt);

  return lik;
}


int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Particle filter DREAM", ' ', "0.1");
  TCLAP::ValueArg<string> json_fn("J","json","JSON input file",true,"","string",cmd);

  TCLAP::UnlabeledMultiArg<string> multi("files","Trees",true,"Input Tree files",false);
  cmd.add(multi);

  try {
    cmd.parse(argc,argv);
  } catch (TCLAP::ArgException &e) { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

  vector<string> fileNames = multi.getValue();

  ifstream in;
  string json_input;
  if (json_fn.getValue() != "") {
    in.open(json_fn.getValue().c_str());
    in.seekg(0,ios::end);
    json_input.reserve(in.tellg());
    in.seekg(0,ios::beg);
    json_input.assign(istreambuf_iterator<char>(in), 
                      istreambuf_iterator<char>());
  }

  // =========================================================================
  // Reading parameters from JSON

  rapidjson::Document jpars;
  try {
    jpars.Parse<0>(json_input.c_str());
    if (! jpars.IsObject()) throw 10;
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  }

  // =========================================================================
  // Reading trees

  Tree tree(fileNames.front().c_str());
  tree.reverse();

  // =========================================================================
  // read DREAM parameters

  dream_pars p;
  dream_pars_default(&p);

  cerr << "Reading parameters..." << flush;
  try {
    dream_pars_read_json(&p,jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }
  cerr << "done" << endl;

  PFPars pf_pars;
  try {
    pf_pars_read_json(&pf_pars,jpars);
  } catch (const char* str) {
    cerr << "PFPars exception: " << str << endl;
  }

  SEISModel::EpiPars seis_pars;
  seis_pars.N     = (p.varHi[0]+p.varLo[0])/2.0;
  seis_pars.beta  = (p.varHi[1]+p.varLo[1])/2.0;
  seis_pars.mu    = (p.varHi[2]+p.varLo[2])/2.0;
  seis_pars.psi   = (p.varHi[3]+p.varLo[3])/2.0;
  seis_pars.gamma = (p.varHi[4]+p.varLo[4])/2.0;
  seis_pars.rho   = 0.0;

  Model* mpt = new SEIS(&seis_pars);
  EpiState* es = new EpiState(SEISModel::nstates);
  (*es)[0] = ((int) seis_pars.N)-1;
  (*es)[1] = 1;
  (*es)[2] = 0;
  (*es)[3] = 1;
  (*es)[4] = 0;
  es->init_branches(tree.max_id()+1);
  es->branches.wake(0);
  es->branches.setCol(0,0);
  // cout << es->to_json() << endl;

  int max_threads = omp_get_max_threads();
  rng::Rng* rng = new rng::GSLRng;
  rng->set_seed(time(NULL));
  rng->alloc(max_threads);

  pf_pars_t lik_pars;
  lik_pars.mpt = mpt;
  lik_pars.tree = &tree;
  lik_pars.pars = &pf_pars;
  lik_pars.rng = rng;
  lik_pars.mpar = &p;
  lik_pars.es = es;

  p.fun = &pf_lik;
  p.funPars = &lik_pars;

  vector<double> last(p.nvar);
  vector<double> proposal(p.nvar);
  vector<double> sigma(p.nvar,0.1);
  last.assign(p.varInit,p.varInit+p.nvar);
  proposal.assign(p.varInit,p.varInit+p.nvar);

  cout << "ll = " << flush;
  double last_lik = pf_lik(proposal.data(),&lik_pars);
  cout << last_lik << endl;

  int accept;
  double new_lik;
  double alpha;
  double drand;
  if (! gsl_isinf(last_lik)) {
    for (int it = 1; it < p.maxEvals; ++it) {
      accept = 0;
      for (size_t j = 0; j < p.nvar; ++j) {
        if (! p.varLock[j]) {
          do {
            (*rng)[0]->gaussian(1,&proposal[j],last[j],sqrt(sigma[j]));
          } while (proposal[j] < p.varLo[j] || proposal[j] > p.varHi[j]);
        }
        cout << setw(12) << last[j] << " ";
      }
      new_lik = pf_lik(proposal.data(),&lik_pars);
      alpha = exp(new_lik-last_lik);
      (*rng)[0]->uniform(1,&drand);
      if (alpha > drand) accept = 1;
      if (accept) {
        last = proposal;
        last_lik = new_lik;
      }
      cout << last_lik << endl;
    }
  }

  dream_pars_free_vars(&p);
  return 0;
}

