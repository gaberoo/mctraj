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

#include "models/Models.h"
using namespace MCTraj;

#include <tclap/CmdLine.h>

#include "Parameters.h"
#include "JSON.h"
#include "pf_pars.h"

int main(int argc, char** argv) {
  // =========================================================================
  // read command line arguments

  TCLAP::CmdLine cmd("Particle filter approximation for marginal tree likelihood", ' ', "0.9.1");
  TCLAP::ValueArg<string> json_fn("J","json","JSON input file",true,"","string",cmd);
  TCLAP::ValueArg<string> out_fn("o","out","output file",false,"","string",cmd);
  TCLAP::ValueArg<string> hist_fn("H","hist","history file",false,"","string",cmd);
  TCLAP::ValueArg<string> _branchfn("B","branchFn","Filename for branch colors",false,"","string",cmd);
  TCLAP::ValueArg<string> _trajfn("C","trajFn","Filename for trajectory",false,"","string",cmd);
  TCLAP::MultiSwitchArg vflag("v","verbose","Increase verbosity",cmd);

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

  try {
    pf_pars_read_json(&pf_pars,jpars);
  } catch (const char* str) {
    cerr << "PFPars exception: " << str << endl;
  }

  pf_pars.vflag += vflag.getValue();

  try {
    p.from_json(jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }

  // =========================================================================

  /*
  TCLAP::ValueArg<int> skip("x","skip","Skip lines of times files",false,0,"string",cmd);
  TCLAP::ValueArg<double> filterTime("f","filter","Min time between filters",false,0.0,"double",cmd);

  TCLAP::SwitchArg _printTraj("O","printTraj","Output trajectory",cmd,false);
  TCLAP::SwitchArg printParticles("P","printParticles","Output particles",cmd,false);
  TCLAP::SwitchArg _fullTree("F","fullTree","Full tree",cmd,false);
  TCLAP::SwitchArg adjZero("z","adjZero","Adjust zero-weight trajectory",cmd,false);
  TCLAP::SwitchArg history("H","history","Store trajectories",cmd,false);
  */

  Tree tree(fileNames.front().c_str(),p.nroot);

  if (tree.getExtant() < 0) {
    cerr << "Cannot run with negative extant lineages at the root." << endl;
    cerr << "0: " << tree.countTypes(0) << endl;
    cerr << "1: " << tree.countTypes(1) << endl;
    return 0;
  }

  // =========================================================================

  size_t npars = p.shifts.size()+1;
  if (npars > 1) {
    if (pf_pars.vflag > 0) cerr << "Rate shifts supplied:" << endl;
    for (size_t i = 0; i < p.shifts.size(); ++i) {
      if (pf_pars.vflag > 0) cerr << "  " << p.shifts[i] << endl;
      tree.addRateShift(p.shifts[i]);
    }
  }

  vector<Pars*> vpars;

  string branch_file = _branchfn.getValue();
  string traj_file = _trajfn.getValue();

  // =========================================================================
  
  if (pf_pars.vflag > 0) {
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
  choose_model(mpt,es,pf_pars,vpars,p);

  double lik = 0.0;

  tree.reverse();

  es->init_branches(tree.max_id()+1,2);
  es->branches.wake(0);
  es->branches.setCol(0,0);

  es->branches.add(0,0);  // overall counter
  es->branches.add(0,1);  // color counter

  // cout << es->to_json() << endl;

  Trajectory* traj = NULL;
  if (pf_pars.print_traj) traj = new Trajectory(*es,mpt);

  ofstream* branch_out = NULL;
  ofstream* traj_out = NULL;

  if (branch_file != "") branch_out = new ofstream(branch_file.c_str());
  if (traj_file != "") traj_out = new ofstream(traj_file.c_str());

  try {
    for (int r = 0; r < pf_pars.reps; ++r) {
      lik = pfLik(mpt,*es,tree,pf_pars,rng,traj);
      cout << lik << endl;
      if (branch_out != NULL) traj->printBranches(*branch_out) << endl;
      if (traj_out != NULL) traj->printFromFirst(*traj_out) << endl;
    }
  } catch (std::exception& e) {
    cout << "EXCEPTION: " << e.what() << endl;
  }

  if (branch_out != NULL) { branch_out->close(); delete branch_out; }
  if (traj_out != NULL) { traj_out->close(); delete traj_out; }

  // if (pf_pars.print_traj) cout << *traj << endl;

  while (vpars.size() > 0) {
    delete vpars.back();
    vpars.pop_back();
  }

  delete traj;
  delete mpt;
  delete es;

  return 0;
}
