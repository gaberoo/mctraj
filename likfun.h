#ifndef __LIKFUN_H__
#define __LIKFUN_H__

#include "MCTraj.h"
#include "pf_pars.h"
#include "pfLik.h"
#include "models/Models.h"
using namespace MCTraj;

#include <cpso/Point.h>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

double pf_likfun(int chain, int gen, const double* state, const void* pars);

inline double pf_likfun_pso(const PSO::Point& state, const void* pars) {
  return pf_likfun(-1,-1,state.data(),pars);
}

#endif
