#include "../pf_pars.h"
#include "../pfLik.h"
#include "../models/Models.h"
using namespace MCTraj;

double pf_i(const PSO::Point& state, const void* pars) 
{
  pf_pars_t& p = *(pf_pars_t*) pars;

  const Pars* oldPt = p.mpt->getPars();
  const IModel::EpiPars* oldPars;

  try {
    oldPars = dynamic_cast<const IModel::EpiPars*>(oldPt);
  } catch (exception& e) {
    cerr << "Error casting pointer!" << endl;
    abort();
  }

  IModel::EpiPars epi(*oldPars);
  epi.beta  = (p.mpar->scale[0] == 'l') ? exp(state[0]) : state[0];
  epi.mu    = (p.mpar->scale[1] == 'l') ? exp(state[1]) : state[1];
  epi.psi   = (p.mpar->scale[2] == 'l') ? exp(state[2]) : state[2];

  if (p.pars->vflag) cerr << epi.to_json() << endl;

  EpiState init(*p.es);

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
//  if (p.otraj != NULL) traj->printFromFirst(*p.otraj) << endl;

  p.mpt->setPars(oldPt);

  return lik;
}

