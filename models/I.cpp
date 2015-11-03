#include "I.h"

double MCTraj::IModel::infRateFun(const EpiState& es, const void* pars, double& trueRate) 
{
  // cerr << "I:INF" << endl;
  EpiPars ep = *(EpiPars*) pars;
  trueRate = ep.beta*es[0];
  if (trueRate < 0.0) trueRate = 0.0;
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::IModel::treeProbInf(const EpiState& es, const void* pars) 
{
  /* Three possibilities upon infection:
   *   (i)   infection involving two tree branches
   *   (ii)  infection involving one tree branch
   *   (iii) infeciton involving zero tree branches
   * The probability that it doesn't involve two tree branches is,
   *   p = 1 - (k choose 2) / (I choose 2) 
   */
  double I = es[0];
  double k = es[1];
  double p = 1.0 - k*(k-1.)/(I*(I-1.));
  return (I >= k) ? p : 0.0;
}

//----------------------------------------------------------------------------

double MCTraj::IModel::recovRateFun(const EpiState& es, const void* pars, double& trueRate) 
{
  // cerr << "I:RECOV" << endl;
  EpiPars ep = *(EpiPars*) pars;
  trueRate = (ep.mu+ep.psi)*es[0];
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::IModel::treeProbRecov(const EpiState& es, const void* pars) 
{
  /* Four possibilities:
   *   (i)   recovery in a non-tree branch
   *   (ii)  sampling in a non-tree branch
   *   (iii) recovery in a tree branch
   *   (iv)  sampling in a tree branch 
   * The probability that it's a recovery in a non-tree branch is
   *   p = (1-s) 
   */
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.mu);
  return (es[0]+1 > es[1]) ? (1.-s) : 0.0;
}

//----------------------------------------------------------------------------

double MCTraj::IModel::treeObsInf(const EpiState& es, const void* pars, double& trueRate) 
{
  EpiPars ep = *(EpiPars*) pars;
  double I = es[0];
  double lambda = ep.beta;
  trueRate = 2.0*lambda/(I+1.);
  return trueRate;
}

//----------------------------------------------------------------------------

double MCTraj::IModel::treeObsRecov(const EpiState& es, const void* pars, double& trueRate) 
{
  EpiPars ep = *(EpiPars*) pars;
  trueRate = (ep.psi > 0) ? ep.psi*es[0] : 1.0;
  return trueRate;
}

/************************************************************************/

double MCTraj::I::sample_rho(const EpiState& es, rng::RngStream* rng, void* pars) const {
  int k = es[1];
  int I = es[0];
  double w = 0.0;
  if (rho() >= 1.0 && I == k) {
    w = 1.0;
  } else if (rho() < 1.0 && I >= k) {
    w = gsl_ran_binomial_pdf(k,rho(),I);
  }
  return w;
}

/************************************************************************/

void MCTraj::I::toTree(const Trajectory& traj,
                       rng::RngStream* rng, 
                       vector<TreeNode>& tree) const 
{
  EpiState x(traj.initState());
  size_t ntrans = traj.transitionCount();
  double time = 0.0;

  const Model* model = traj.getModel();
  vector< vector<int> > states(model->n());

  size_t i, j;
  int branch_id = 0;
  int epi_id = 0;

  int lineageStates[] = { 1, 0 };

  // Initialize
  for (i = 0; i < x.numStates(); ++i) {
    if (lineageStates[i]) {
      for (j = 0; j < (size_t) x[i]; ++j) {
        tree.push_back(TreeNode());
        branch_id = tree.size()-1;
        tree.back().code = branch_id;
        tree.back().epi_id = epi_id++;
        tree.back().epi_state = i;
        states[i].push_back(branch_id);
      }
    }
  }

  int rand;
  int parent;

  double p = 1.0;
  double r;
  StateTransition st;

  for (size_t i = 0; i < ntrans; ++i) {
    st = traj.getTrans(i); 
    x += st.getTrans();
    time += st.atTime();

    switch (st.etype()) {
      case 1:
        rng->uniform_int(1,&rand,0,states[0].size());

        parent = states[0][rand];
        states[0].erase(states[0].begin()+rand);

        tree[parent].age = time;
        tree[parent].n_off = 2;

        // old infected
        tree.push_back(TreeNode());
        branch_id = tree.size()-1;

        tree[branch_id].code = branch_id;
        tree[branch_id].parent = tree[parent].code;
        tree[branch_id].age = -1.0;
        tree[branch_id].ancestor_age = tree[parent].age;
        tree[branch_id].epi_id = tree[parent].epi_id;
        tree[branch_id].epi_state = tree[parent].epi_state;

        states[0].push_back(branch_id);
        tree[parent].off.push_back(branch_id);

        // new infected
        tree.push_back(TreeNode());
        branch_id = tree.size()-1;

        tree[branch_id].code = branch_id;
        tree[branch_id].parent = tree[parent].code;
        tree[branch_id].age = -1.0;
        tree[branch_id].ancestor_age = tree[parent].age;
        tree[branch_id].epi_id = epi_id++;
        tree[branch_id].epi_state = 1;

        states[0].push_back(branch_id);
        tree[parent].off.push_back(branch_id);

        break;

      case 0:
        rng->uniform_int(1,&rand,0,states[0].size());
        parent = states[0][rand];
        states[0].erase(states[0].begin()+rand);
        tree[parent].age = time;
        tree[parent].n_off = 0;

        p = st.getType()->applyProb(x,model->getPars());
        rng->uniform(1,&r);

        if (r < p) {
          add_extant(tree,parent);
          // tree[parent].extant_off = 2;
        }

        break;
    }
  }

  // print the remaining branches
  for (i = 0; i < x.numStates(); ++i) {
    if (lineageStates[i]) {
      for (j = 0; j < states[i].size(); ++j) {
        tree[states[i][j]].age = time;
      }
    }
  }

}

