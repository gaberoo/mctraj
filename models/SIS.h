#ifndef __SISMODEL_H__
#define __SISMODEL_H__

#include <gsl/gsl_math.h>
#include "../Model.h"
#include "../EpiState.h"

namespace MCTraj {
  namespace SISModel {
    struct EpiPars {
      double N;
      double beta;
      double mu;
      double psi;
      double rho;
    };

    /* S, I(all), I(tree) */
    const size_t nstates = 3;

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { 1, -1, 0 };

    /************************************************************************/

    double treeObsInf(const EpiState& es, const void* pars);
    const int obsInfChange[] = { -1, 1, 1 };

    /************************************************************************/

    double treeObsRecov(const EpiState& es, const void* pars);
    const int obsRecovChange[] = { 1, -1, -1 };
  }

  /**************************************************************************/

  class SIS : public Model {
    public:
      SIS(const SISModel::EpiPars* p) {
        nstates = SISModel::nstates;
        pars = p;
        rho = p->rho;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov",
                                                SISModel::nstates,
                                                SISModel::recoverChange,
                                                SISModel::recovRateFun,
                                                SISModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov",
                                              SISModel::nstates,
                                              SISModel::obsRecovChange,
                                              SISModel::treeObsRecov,
                                              oneProb));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf",
                                                SISModel::nstates,
                                                SISModel::infChange,
                                                SISModel::infRateFun,
                                                SISModel::treeProbInf));

        obsTypes.push_back(new TransitionType("obsInf",
                                              SISModel::nstates,
                                              SISModel::obsInfChange,
                                              SISModel::treeObsInf,
                                              oneProb));
        simEvent.push_back(1);
      }

      virtual ~SIS() {}

      double sample_rho(const EpiState& es, gsl_rng* rng, void* pars = NULL) const;
  };
}

#endif // __SISMODEL_H__
