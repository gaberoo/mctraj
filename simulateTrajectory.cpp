#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "MCTraj.h"
#include "Model.h"
#include "models/Models.h"
using namespace MCTraj;

#ifdef MKLRNG
#include <rng/MKLRng.h>
#else
#include <rng/GSLRng.h>
#endif

#include <tclap/CmdLine.h>
#include "JSON.h"

int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Simulate trees using particle filter models", ' ', "0.9.1");

  TCLAP::ValueArg<string> json_fn("J","json","JSON input file",true,"","string",cmd);
  TCLAP::MultiSwitchArg vflag("v","verbose","Increase verbosity",cmd);
  TCLAP::ValueArg<double> mT("T","maxTime","Maximum simulation time",false,100.0,"double",cmd);

  try {
    cmd.parse(argc,argv);
  } catch (TCLAP::ArgException &e) { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

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

  // Setup random number generators
  rng::RngStream* rng = NULL;
#ifdef MKLRNG
  rng = new rng::MKLStream;
#else
  rng = new rng::GSLStream;
#endif
  rng->alloc(pf_pars.seed);

  size_t npars = p.shifts.size()+1;
  if (npars > 1) {
    if (pf_pars.vflag > 0) cerr << "Rate shifts supplied:" << endl;
    for (size_t i = 0; i < p.shifts.size(); ++i) {
      if (pf_pars.vflag > 0) cerr << "  " << p.shifts[i] << endl;
    }
  }
  vector<Pars*> vpars;

  EpiState* es = NULL;
  Model* mpt = NULL;
  choose_model(mpt,es,pf_pars,vpars,p);

  Trajectory traj(*es,mpt);

  p.shifts.push_back(mT.getValue());
  sort(p.shifts.begin(),p.shifts.end());
  size_t ti = 0;
  do {
    traj.simulateTrajectory(p.shifts[ti],&vpars[ti],rng,false);
  } while (++ti < p.shifts.size());

  vector<TreeNode> tree;
  vector<TreeNode> phylo;

  mpt->toTree(traj,rng,tree);

  // for (size_t i = 0; i < tree.size(); ++i) cout << tree[i] << endl;

//  string newick = to_newick(tree,0);
//  string sample_newick;
//  only_sampled(tree,sample_newick,0);
//  cout << "begin trees;" << endl;
//  cout << "tree 'full_tree' = " << newick << ";" << endl;
//  cout << "tree 'sampled_tree' = " << sample_newick << ";" << endl;
//  cout << "end;" << endl;

  only_sampled(tree,phylo,0);

  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
  json_w.StartObject(); {
    json_w.String("model"); vpars[0]->json(json_w);
    json_w.String("traj"); traj.json(json_w);
    json_w.String("full_tree"); json_w.StartArray(); {
      for (size_t i = 0; i < tree.size(); ++i) tree[i].json(json_w);
    } json_w.EndArray();
    json_w.String("sampled"); json_w.StartArray(); {
      for (size_t i = 0; i < phylo.size(); ++i) phylo[i].json(json_w);
    } json_w.EndArray();
  } json_w.EndObject();

  cout << buf.GetString() << endl;

  while (vpars.size() > 0) { delete vpars.back(); vpars.pop_back(); }
  delete mpt;
  delete es;

  return 0;
}
