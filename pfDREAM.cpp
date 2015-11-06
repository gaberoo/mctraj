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

#include "JSON.h"
#include "pf_pars.h"
#include "pso_pars.h"

#include "likfun.h"

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

  string json_input = read_json(json_fn.getValue());

  rapidjson::Document jpars;
  get_json(jpars,json_input);

  PFPars pf_pars;
  Parameters pars;

  try {
    pf_pars.from_json(jpars);
  } catch (const char* str) {
    cerr << "PFPars exception: " << str << endl;
  }

  try {
    pars.from_json(jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
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

#if defined(_OPENMP)
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

  rng::Rng* rng = new rng::GSLRng;
  rng->set_seed(time(NULL));
  rng->alloc(max_threads);

  int nroot = pars.nroot;

  IModel::EpiPars       i_pars;
  SISModel::EpiPars   sis_pars;
  SEISModel::EpiPars seis_pars;

  Tree tree(fileNames.front().c_str(),nroot);
  tree.reverse();

  // =========================================================================

  vector<Pars*> vpars;
  EpiState* es = NULL;
  Model* mpt = NULL;
  choose_model(mpt,es,pf_pars,vpars,pars);

  // =========================================================================

  string pf_hist_fn = "";
  rapidjson::Value::MemberIterator it = jpars.FindMember("pf");
  if (it != jpars.MemberEnd()) {
    rapidjson::Value& d = it->value;
    it = d.FindMember("hist_fn");
    if (it != d.MemberEnd()) pf_hist_fn = it->value.GetString();
  }

  pf_pars_t lik_pars;

  lik_pars.mpt = mpt;
  lik_pars.tree = &tree;
  lik_pars.pars = &pf_pars;
  lik_pars.rng = rng;
  lik_pars.pfpar = &pars;
  lik_pars.dpar = &p;
  lik_pars.es = es;

  if (pf_hist_fn != "") lik_pars.otraj = new ofstream(pf_hist_fn.c_str());

  p.fun = &pf_likfun;
  p.funPars = &lik_pars;

//  double state[] = { 100.0, 0.0, 0.1, 0.4, -2.3 };
//  cout << pf_pars.vflag << " > ll = " << flush;
//  double ll = pf_lik(state,&lik_pars);
//  cout << ll << endl;

  dream(&p,(*rng)[0]);

  dream_pars_free_vars(&p);

  return 0;
}

