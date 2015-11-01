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

#include "pf_pars.h"
#include "pso_pars.h"

#include "funcs/I.h"
#include "funcs/SIS.h"

int main(int argc, char** argv) {
  // =========================================================================
  // read command line arguments

  TCLAP::CmdLine cmd("Particle filter DREAM", ' ', "0.2");
  TCLAP::ValueArg<string> json_fn("J","json","JSON input file",true,"","string",cmd);

  TCLAP::UnlabeledMultiArg<string> multi("files","Trees",true,"Input Tree files",false);
  cmd.add(multi);

  try {
    cmd.parse(argc,argv);
  } catch (TCLAP::ArgException &e) { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

  vector<string> fileNames = multi.getValue();

  // =========================================================================
  // read JSON file

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

#if defined(_OPENMP)
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

  rng::Rng* rng = new rng::GSLRng;
  rng->set_seed(time(NULL));
  rng->alloc(max_threads);

  int nroot = 0;
  it = jpars.FindMember("nroot");
  if (it != jpars.MemberEnd()) nroot = it->value.GetInt();

  IModel::EpiPars       i_pars;
  SISModel::EpiPars   sis_pars;
  SEISModel::EpiPars seis_pars;

  Tree tree(fileNames.front().c_str(),nroot);
  tree.reverse();

  // =========================================================================

  EpiState* es = NULL;
  Model* mpt = NULL;

  int kase = 0;

  switch (pf_pars.model_type) {
    case 'I':
    case 0:
      if (ppso.vflag > 1) cerr << "I model." << endl;
      i_pars.beta = p.initVar(0);
      i_pars.mu   = p.initVar(1);
      i_pars.psi  = p.initVar(2);
      i_pars.rho  = p.initVar(3);
      mpt = new I(&i_pars);
      es = new EpiState(IModel::nstates);
      (*es)[0] = nroot;
      (*es)[1] = nroot;

      p.fun = &pf_i;

      break;

    case 'R':
    case 2:
      kase = 1;
    case 'S':
    case 1:
    default:
      if (ppso.vflag > 1) cerr << "SIS/SIR model." << endl;
      sis_pars.N    = p.initVar(0);
      sis_pars.beta = p.initVar(1);
      sis_pars.mu   = p.initVar(2);
      sis_pars.psi  = p.initVar(3);
      sis_pars.rho  = p.initVar(4);

      if (kase == 1) mpt = new SIR(&sis_pars);
      else mpt = new SIS(&sis_pars);

      es = new EpiState(SISModel::nstates);
      (*es)[0] = (int) sis_pars.N - nroot;
      (*es)[1] = nroot;
      (*es)[2] = nroot;

      p.fun = &pf_sis;

      break;

    case 'E':
    case 3:
      cerr << "SEIS not yet implemented!" << endl;
      break;
  }

  pf_pars_t lik_pars;
  lik_pars.mpt = mpt;
  lik_pars.tree = &tree;
  lik_pars.pars = &pf_pars;
  lik_pars.rng = rng;
  lik_pars.mpar = &p;
  lik_pars.es = es;
  if (pf_hist_fn != "") lik_pars.otraj = new ofstream(pf_hist_fn.c_str());

  p.funPars = &lik_pars;

//  double state[] = { 100.0, 0.0, 0.1, 0.4, -2.3 };
//  cout << pf_pars.vflag << " > ll = " << flush;
//  double ll = pf_lik(state,&lik_pars);
//  cout << ll << endl;

  dream(&p,(*rng)[0]);

  dream_pars_free_vars(&p);
  return 0;
}

