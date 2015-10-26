#include "../pf_pars.h"
#include "../pfLik.h"
#include "../models/Models.h"
using namespace MCTraj;

double pf_seis(const PSO::Point& state, const void* pars) 
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

  if (p.pars->vflag > 1) cerr << "Starting calculation." << endl;
  double lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,traj);

  if (p.pars->vflag > 1) cerr << "Finished calculation." << endl;

//  if (p.obranch != NULL) traj->printBranches(*p.obranch) << endl;
//  if (p.otraj != NULL) traj->printFromFirst(*p.otraj) << endl;

  p.mpt->setPars(oldPt);

  return lik;
}


