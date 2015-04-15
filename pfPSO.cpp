#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <tclap/CmdLine.h>
#include <rapidjson/document.h>

#include <rng/GSLRng.h>
#include <cpso/Swarm.h>

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
  const PSO::Parameters* mpar;
  ostream* otraj;
  ostream* obranch;
} pf_pars_t;

double pf_lik(const PSO::Point& state, const void* pars) 
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

  Trajectory* traj = NULL;
  if (p.otraj != NULL) traj = new Trajectory(init,p.mpt);

  double lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,traj);

  if (p.obranch != NULL) traj->printBranches(*p.obranch) << endl;
  if (p.otraj != NULL) traj->printFromFirst(*p.otraj) << endl;

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

  PSO::Parameters p;

  cerr << "Reading parameters..." << flush;
  try {
    p.from_json(jpars);
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
  seis_pars.N     = (p.ub[0]+p.lb[0])/2.0;
  seis_pars.beta  = (p.ub[1]+p.lb[1])/2.0;
  seis_pars.mu    = (p.ub[2]+p.lb[2])/2.0;
  seis_pars.psi   = (p.ub[3]+p.lb[3])/2.0;
  seis_pars.gamma = (p.ub[4]+p.lb[4])/2.0;
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

  p.evalFunc = &pf_lik;
  p.evalParams = &lik_pars;

  // =========================================================================
  // PSO parameters

  int max_evals = 100;
  int num_particles = 20;
  int init_type = 1;
  int vflag = 0;
  int slowdown = 0;
  ostream* out = NULL;
  ostream* hist = NULL;

  rapidjson::Value::MemberIterator it = jpars.FindMember("pso");
  if (it != jpars.MemberEnd()) {
    rapidjson::Value& d = it->value;
    it = d.FindMember("max_evals");
    if (it != d.MemberEnd()) max_evals = it->value.GetInt();
    it = d.FindMember("num_particles");
    if (it != d.MemberEnd()) num_particles = it->value.GetInt();
    it = d.FindMember("init_type");
    if (it != d.MemberEnd()) init_type = it->value.GetInt();
    it = d.FindMember("vflag");
    if (it != d.MemberEnd()) vflag = it->value.GetInt();
    it = d.FindMember("slowdown");
    if (it != d.MemberEnd()) slowdown = it->value.GetInt();
    it = d.FindMember("out_fn");
    if (it != d.MemberEnd()) out = new ofstream(it->value.GetString());
    it = d.FindMember("hist_fn");
    if (it != d.MemberEnd()) hist = new ofstream(it->value.GetString());
  }

  // =========================================================================
  // Run PSO

  PSO::Swarm s(num_particles,5,&p);

  PSO::Point phi_p(5,1.49445);
  PSO::Point phi_g(5,1.49445);
  PSO::Point omega(5,0.729);
  s.setVars(phi_p,phi_g,omega);

  s.begin(argc,argv);
  s.initialize(init_type);
  s.evaluate();
  s.display();
#ifdef USE_MPI
  s.run_mpi(max_evals,vflag,out,hist);
#else
  s.run(max_evals,slowdown,vflag,out,hist);
#endif
  s.end();

  return 0;
}

