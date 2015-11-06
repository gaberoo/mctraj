#ifndef __PF_PARS_H__
#define __PF_PARS_H__

#include "Parameters.h"
#include <cpso/Parameters.h>
#include <cdream/dream.h>

typedef struct {
  Model* mpt = NULL;
  const Tree* tree = NULL;
  const PFPars* pars = NULL;
  const EpiState* es = NULL;
  rng::Rng* rng = NULL;
  const Parameters* pfpar = NULL;
  const PSO::Parameters* mpar = NULL;
  const dream_pars* dpar = NULL;
  ostream* otraj = NULL;
  ostream* obranch = NULL;
} pf_pars_t;

// ===========================================================================

#endif

