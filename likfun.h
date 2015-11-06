#ifndef __LIKFUN_H__
#define __LIKFUN_H__

#include "MCTraj.h"
#include "pf_pars.h"
#include "pfLik.h"
#include "models/Models.h"
using namespace MCTraj;

#include <cpso/Point.h>

double pf_likfun(const double* state, const void* pars);

inline double pf_likfun_pso(const PSO::Point& state, const void* pars) {
  return pf_likfun(state.data(),pars);
}

#endif
