#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "MCTraj.h"
#include "Model.h"
#include "models/SEIS.h"
using namespace MCTraj;

#ifdef MKLRNG
#include <rng/MKLRng.h>
#else
#include <rng/GSLRng.h>
#endif

int main(int argc, char** argv) {
  SEISModel::EpiPars seis_pars;
  seis_pars.N     = 100.0;
  seis_pars.beta  = 1.0;
  seis_pars.mu    = 0.1;
  seis_pars.psi   = 0.4;
  seis_pars.rho   = 0.0;
  seis_pars.gamma = 0.1;
  double maxTime = 100.0;

  SEIS model(&seis_pars);

  EpiState init(SEISModel::nstates);
  init[0] = seis_pars.N-1;
  init[1] = 0;
  init[2] = 1;
  init[3] = 0;
  init[4] = 0;

  Trajectory traj(init,&model);

  unsigned long seed = (argc > 1) ? atoi(argv[1]) : time(NULL);

  // Setup random number generators
  rng::RngStream* rng = NULL;
#ifdef MKLRNG
  rng = new rng::MKLStream;
#else
  rng = new rng::GSLStream;
#endif
  rng->alloc(seed);

  traj.simulateTrajectory(maxTime,&seis_pars,rng);

  vector<TreeNode> tree;
  vector<TreeNode> phylo;

  int lineageStates[] = { 0, 1, 1, 0, 0 };
  traj.toTree(rng,tree,lineageStates);

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
    json_w.String("model"); seis_pars.json(json_w);
    json_w.String("full_tree"); json_w.StartArray(); {
      for (size_t i = 0; i < tree.size(); ++i) tree[i].json(json_w);
    } json_w.EndArray();
    json_w.String("sampled"); json_w.StartArray(); {
      for (size_t i = 0; i < phylo.size(); ++i) phylo[i].json(json_w);
    } json_w.EndArray();
  } json_w.EndObject();

  cout << buf.GetString() << endl;

  return 0;
}
