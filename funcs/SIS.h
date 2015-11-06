#include "../pf_pars.h"
#include "../pfLik.h"
#include "../models/Models.h"
using namespace MCTraj;

double pf_sis(const double* state, const void* pars);

inline double pf_sis_pso(const PSO::Point& state, const void* pars) {
  return pf_sis(state.data(),pars);
}

double pf_sis(const double* state, const void* pars) 
{
  pf_pars_t& p = *(pf_pars_t*) pars;

  const Pars* oldPt = p.mpt->getPars();
  const SISModel::EpiPars* oldPars;

  try {
    oldPars = dynamic_cast<const SISModel::EpiPars*>(oldPt);
  } catch (exception& e) {
    cerr << "Error casting pointer!" << endl;
    abort();
  }

  const char* scale = NULL;
  if (p.mpar != NULL) {
    scale = p.mpar->scale.data();
  } else if (p.dpar != NULL) {
    scale = p.dpar->scale;
  } else {
    return -INFINITY;
  }

  SISModel::EpiPars epi(*oldPars);
  epi.N     = (scale[0] == 'l') ? exp(state[0]) : state[0];
  epi.beta  = (scale[1] == 'l') ? exp(state[1]) : state[1];
  epi.mu    = (scale[2] == 'l') ? exp(state[2]) : state[2];
  epi.psi   = (scale[3] == 'l') ? exp(state[3]) : state[3];

  if (p.pars->vflag) cerr << epi.to_json() << endl;

  EpiState init(*p.es);
  init[0] = (int) epi.N - init[1]; 
  if (init[0] < 0) return -INFINITY;

  /*
  init.init_branches(p.tree->max_id()+1,2);
  init.branches.wake(0);
  init.branches.setCol(0,0);

  init.branches.add(0,0);  // overall counter
  init.branches.add(0,1);  // color counter
  */

  Trajectory* traj = NULL;
  if (p.otraj != NULL) traj = new Trajectory(init,p.mpt);

  p.mpt->setPars(&epi);

  double lik = -INFINITY;

  try {
    if (p.pars->vflag > 1) cerr << "Starting calculation." << endl;
    lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,traj);
    if (p.pars->vflag > 1) cerr << "Finished calculation." << endl;
  } catch (std::exception& e) {
    cout << "EXCEPTION: " << e.what() << endl;
  }

//  if (p.obranch != NULL) traj->printBranches(*p.obranch) << endl;
  if (p.otraj != NULL) *p.otraj << traj->to_json() << endl;

  p.mpt->setPars(oldPt);

  return lik;
}


