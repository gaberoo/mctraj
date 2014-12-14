#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "Model.h"
#include "models/SEIS.h"
using namespace MCTraj;

int main(int argc, char** argv) {
  SEISModel::EpiPars seis_pars;
  seis_pars.N     = 100;
  seis_pars.beta  = 1.0;
  seis_pars.mu    = 0.1;
  seis_pars.psi   = 0.1;
  seis_pars.rho   = 0.0;
  seis_pars.gamma = 0.01;

  double maxTime = 100.0;

  SEIS model(&seis_pars);

  EpiState init(SEISModel::nstates);
  init[0] = 99;
  init[1] = 0;
  init[2] = 1;
  init[3] = 0;
  init[4] = 0;

  Trajectory traj(init,&model);

  unsigned long seed = (argc > 1) ? atoi(argv[1]) : time(NULL);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,seed);

  traj.simulateTrajectory(maxTime,&seis_pars,rng);

  vector<TreeNode> tree;
  vector<TreeNode> phylo;
  traj.toTree(rng,tree);

  // for (size_t i = 0; i < tree.size(); ++i) cout << tree[i] << endl;

//  string newick = to_newick(tree,0);
//  string sample_newick;
//  only_sampled(tree,sample_newick,0);
//  cout << "begin trees;" << endl;
//  cout << "tree 'full_tree' = " << newick << ";" << endl;
//  cout << "tree 'sampled_tree' = " << sample_newick << ";" << endl;
//  cout << "end;" << endl;

  only_sampled(tree,phylo,0);
  cout << "[";
  for (size_t i = 0; i < phylo.size(); ++i) {
    cout << phylo[i];
    if (i < phylo.size()-1) cout << ",";
    cout << endl;
  }
  cout << "]" << endl;

//  cout << setw(5)  << "\"P\"" << " "
//      << setw(5)  << "\"ID\"" << " "
//      << setw(12) << "\"age\"" << " "
//      << setw(5)  << "\"n_off\"" << " "
//      << setw(5)  << "\"exoff\"" << " "
//      << setw(12) << "\"P_age\"" << " "
//      << setw(5)  << "\"eID\"" << " "
//      << setw(5)  << "\"state\"" << endl;
//
//  vector<TreeNode> sample;
//  only_sampled(tree,sample,0);
//  for (size_t i = 0; i < sample.size(); ++i) {
//    cout << sample[i] << endl;
//  }

  gsl_rng_free(rng);

  return 0;
}
