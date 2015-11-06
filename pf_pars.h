#ifndef __PF_PARS_H__
#define __PF_PARS_H__

#include "MCTraj.h"
#include "Tree.h"
#include "Parameters.h"
#include "Trajectory.h"
#include <cpso/Parameters.h>
#include <cdream/dream.h>

namespace MCTraj {
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
    string scales = "";
  } pf_pars_t;
}

// ===========================================================================

#endif

