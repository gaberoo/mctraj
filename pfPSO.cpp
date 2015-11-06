#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <tclap/CmdLine.h>
#include <rapidjson/document.h>

#include <rng/GSLRng.h>
#include <cpso/Swarm.h>

#include "pfLik.h"
#include "models/Models.h"
using namespace MCTraj;

#include "JSON.h"
#include "pf_pars.h"
#include "pso_pars.h"

#include "funcs/I.h"
#include "funcs/SIS.h"

int main(int argc, char** argv) {
  // =========================================================================
  // read command line arguments

  TCLAP::CmdLine cmd("Particle filter PSO", ' ', "0.2");
  TCLAP::ValueArg<string> json_fn("J","json","JSON input file",true,"","string",cmd);
  TCLAP::ValueArg<string> out_fn("o","out","output file",false,"","string",cmd);
  TCLAP::ValueArg<string> hist_fn("H","hist","history file",false,"","string",cmd);

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
  Parameters p;
  pso_pars_t ppso;

  try {
    pf_pars_read_json(&pf_pars,jpars);
  } catch (const char* str) {
    cerr << "PFPars exception: " << str << endl;
  }

  try {
    p.from_json(jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }

  // PSO parameters

  PSO::Parameters pso_pars;
  try {
    pso_pars.from_json(jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }

  rapidjson::Value::MemberIterator it = jpars.FindMember("pso");
  if (it != jpars.MemberEnd()) ppso.from_json(it->value);
  if (hist_fn.getValue() != "") {
    delete ppso.hist;
    ppso.hist = new ofstream(hist_fn.getValue().c_str());
  }
  if (out_fn.getValue() != "") {
    delete ppso.out;
    ppso.out = new ofstream(out_fn.getValue().c_str());
  }

  // =========================================================================

#if defined(_OPENMP)
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

  rng::Rng* rng = new rng::GSLRng;
  rng->set_seed(time(NULL));
  rng->alloc(max_threads);

  int nroot = p.nroot;

  string pf_hist_fn = "";
  it = jpars.FindMember("pf");
  if (it != jpars.MemberEnd()) {
    rapidjson::Value& d = it->value;
    it = d.FindMember("hist_fn");
    if (it != d.MemberEnd()) pf_hist_fn = it->value.GetString();
  }

  IModel::EpiPars       i_pars;
  SISModel::EpiPars   sis_pars;
  SEISModel::EpiPars seis_pars;

  Tree tree(fileNames.front().c_str(),nroot);
  tree.reverse();

  // =========================================================================

  vector<Pars*> vpars;
  EpiState* es = NULL;
  Model* mpt = NULL;
  choose_model(mpt,es,pf_pars,vpars,p);

  switch (pf_pars.model_type) {
    case 'I':
    case 0:
      pso_pars.evalFunc = &pf_i_pso;
      break;

    case 'S':
    case 1:
    case 'R':
    case 2:
    default:
      pso_pars.evalFunc = &pf_sis_pso;
      break;
  }

  pf_pars_t lik_pars;
  lik_pars.mpt  = mpt;
  lik_pars.tree = &tree;
  lik_pars.pars = &pf_pars;
  lik_pars.rng  = rng;
  lik_pars.mpar = &pso_pars;
  lik_pars.es   = es;
  if (pf_hist_fn != "") lik_pars.otraj = new ofstream(pf_hist_fn.c_str());

  pso_pars.evalParams = &lik_pars;

  // =========================================================================
  // Run PSO

  PSO::Swarm s(ppso.num_particles,5,&pso_pars);

  PSO::Point phi_p(5,1.49445);
  PSO::Point phi_g(5,1.49445);
  PSO::Point omega(5,0.729);
  s.setVars(phi_p,phi_g,omega);

  s.begin(argc,argv);
  s.initialize(ppso.init_type);
  s.evaluate();
  s.display();
#ifdef USE_MPI
  s.run_mpi(ppso.max_evals,ppso.vflag,ppso.out,ppso.hist);
#else
  s.run(ppso.max_evals,ppso.slowdown,ppso.vflag,ppso.out,ppso.hist);
#endif
  s.end();

  if (lik_pars.otraj != NULL) delete lik_pars.otraj;

  return 0;
}

