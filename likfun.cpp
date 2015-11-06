#include "likfun.h"

double pf_likfun(const double* state, const void* pars) 
{
  double lik = -INFINITY;

  pf_pars_t& p = *(pf_pars_t*) pars;

  string scales = p.pfpar->scales();
  if (p.mpar != NULL) scales.assign(p.mpar->scale.begin(),p.mpar->scale.end());
  else if (p.dpar != NULL) scales = p.dpar->scale;

  p.mpt->getPars()->from_state(state,scales.c_str());
  if (p.pars->vflag) cerr << p.mpt->getPars()->to_json() << endl;
  EpiState init(p.mpt->initState(*(p.pfpar)));

  Trajectory* traj = NULL;
  if (p.otraj != NULL) traj = new Trajectory(init,p.mpt);

  try {
    if (p.pars->vflag > 1) cerr << "Starting calculation." << endl;
    lik = pfLik(p.mpt,init,*(p.tree),*(p.pars),p.rng,traj);
    if (p.pars->vflag > 1) cerr << "Finished calculation." << endl;
  } catch (std::exception& e) {
    cout << "EXCEPTION: " << e.what() << endl;
  }

//  if (p.obranch != NULL) traj->printBranches(*p.obranch) << endl;
  if (p.otraj != NULL) *p.otraj << traj->to_json() << endl;

  return lik;
}


