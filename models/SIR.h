#ifndef __SIRMODEL_H__
#define __SIRMODEL_H__

#include <rng/Rng.h>
#include "../Model.h"
#include "SIS.h"
#include "../EpiState.h"

namespace MCTraj {
  namespace SIRModel {
    typedef SISModel::EpiPars EpiPars;

    /* S, I(all), I(tree) */
    const size_t nstates = 4;

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0, 0 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { 0, -1, 0, 1 };

    /************************************************************************/

    double treeObsInf(const EpiState& es, const void* pars);
    const int obsInfChange[] = { -1, 1, 1, 0 };

    /************************************************************************/

    double treeObsRecov(const EpiState& es, const void* pars);
    const int obsRecovChange[] = { 0, -1, -1, 1 };
  }

  /**************************************************************************/

  class SIR : public Model {
    public:
      SIR(const SIRModel::EpiPars* p) {
        nstates = SIRModel::nstates;
        pars = p;
        rho = p->rho;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov",
                                                SIRModel::nstates,
                                                SIRModel::recoverChange,
                                                SIRModel::recovRateFun,
                                                SIRModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov",
                                              SIRModel::nstates,
                                              SIRModel::obsRecovChange,
                                              SIRModel::treeObsRecov,
                                              oneProb));
        simEvent.push_back(1);


        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf",
                                                SIRModel::nstates,
                                                SIRModel::infChange,
                                                SIRModel::infRateFun,
                                                SIRModel::treeProbInf));

        obsTypes.push_back(new TransitionType("obsInf",
                                              SIRModel::nstates,
                                              SIRModel::obsInfChange,
                                              SIRModel::treeObsInf,
                                              oneProb));
        simEvent.push_back(1);

      }

      virtual ~SIR() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      inline bool validState(const EpiState& es) const { return true; }
  };
}

#endif // __SIRMODEL_H__
