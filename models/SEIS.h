#ifndef __SEISMODEL_H__
#define __SEISMODEL_H__

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "../Model.h"
#include "../EpiState.h"
#include "../StateTransition.h"

namespace MCTraj {
  namespace SEISModel {
    struct EpiPars {
      double N;
      double beta;
      double mu;
      double psi;
      double rho;
      double gamma;
    };

    const size_t nstates = 5; /* S, E, I, kE, kI */

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0, 0, 0 };

    /************************************************************************/

    double transRateFun(const EpiState& es, const void* pars);
    double treeProbTrans(const EpiState& es, const void* pars);
    const int transChange[] = { 0, -1, 1, 0, 0 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { 1, 0, -1, 0, 0 };

    /************************************************************************/

    double treeObsInf(const EpiState& es, const void* pars);
    const int obsInfChange[] = { -1, 1, 0, 1, 0 };

    /************************************************************************/

    double treeObsTrans(const EpiState& es, const void* pars);
    const int obsTransChange[] = { 0, -1, 1, -1, 1 };

    /************************************************************************/

    double treeObsRecov(const EpiState& es, const void* pars);
    const int obsRecovChange[] = { 1, 0, -1, 0, -1 };

    /************************************************************************/

    int branchInf(const EpiState& es, gsl_rng* rng, 
                   StateTransition& st, const void* pars);
    int branchTrans(const EpiState& es, gsl_rng* rng, 
                     StateTransition& st, const void* pars);

    int obsBranchInf(const EpiState& es, gsl_rng* rng, 
                      StateTransition& st, const void* pars);
    int obsBranchRecov(const EpiState& es, gsl_rng* rng, 
                        StateTransition& st, const void* pars);
    int obsBranchTrans(const EpiState& es, gsl_rng* rng, 
                        StateTransition& st, const void* pars);
  }

  /**************************************************************************/

  class SEIS : public Model {
    public:
      SEIS(const SEISModel::EpiPars* p) {
        nstates = SEISModel::nstates;
        pars = p;
        rho = p->rho;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov",
                                                SEISModel::nstates,
                                                SEISModel::recoverChange,
                                                SEISModel::recovRateFun,
                                                SEISModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov",
                                              SEISModel::nstates,
                                              SEISModel::obsRecovChange,
                                              SEISModel::treeObsRecov,
                                              oneProb,
                                              SEISModel::obsBranchRecov));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf",
                                                SEISModel::nstates,
                                                SEISModel::infChange,
                                                SEISModel::infRateFun,
                                                SEISModel::treeProbInf,
                                                SEISModel::branchInf));

        obsTypes.push_back(new TransitionType("obsInf",
                                              SEISModel::nstates,
                                              SEISModel::obsInfChange,
                                              SEISModel::treeObsInf,
                                              oneProb,
                                              SEISModel::obsBranchInf));
        simEvent.push_back(1);

        /* transition events */
        typeMap[70] = 2;
        transTypes.push_back(new TransitionType("simTrans",
                                                SEISModel::nstates,
                                                SEISModel::transChange,
                                                SEISModel::transRateFun,
                                                SEISModel::treeProbTrans,
                                                SEISModel::branchTrans));

        obsTypes.push_back(new TransitionType("obsTrans",
                                              SEISModel::nstates,
                                              SEISModel::obsTransChange,
                                              SEISModel::treeObsTrans,
                                              oneProb,
                                              SEISModel::obsBranchTrans));
        simEvent.push_back(1);

      }

      virtual ~SEIS() {}
  };
}

#endif // __SEISMODEL_H__
