#ifndef __PF_PARS_H__
#define __PF_PARS_H__

typedef struct {
  Model* mpt;
  const Tree* tree;
  const PFPars* pars;
  const EpiState* es;
  rng::Rng* rng;
  const PSO::Parameters* mpar;
  ostream* otraj;
  ostream* obranch;
} pf_pars_t;

#endif

