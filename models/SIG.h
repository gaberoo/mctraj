#ifndef __SIGMODEL_H__
#define __SIGMODEL_H__

#include <rng/Rng.h>
#include "../Model.h"
#include "SIS.h"
#include "../EpiState.h"

namespace MCTraj {
  namespace SIGModel {
    typedef SISModel::EpiPars EpiPars;

    /* Susceptible, Infected, Sampled, Recovered, Lineage */
    const size_t nstates = 5;

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0, 0, 0 };

    double treeObsInf(const EpiState& es, const void* pars);
    const int obsInfChange[] = { -1, 1, 0, 0, 1 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { 0, -1, 0, 1, 0 };

    double treeObsRecov(const EpiState& es, const void* pars);
    const int obsRecovChange[] = { 0, -1, 0, 1, -1 };

    /************************************************************************/

    double sampRateFun(const EpiState& es, const void* pars);
    double treeProbSample(const EpiState& es, const void* pars);
    const int sampleChange[] = { 0, -1, 1, 0, 0 };

    double treeObsSample(const EpiState& es, const void* pars);
    const int obsSampleChange[] = { 0, -1, 1, 0, -1 };
  }

  /**************************************************************************/

  class SIG : public Model {
    public:
      SIG(const SIG& m) : Model(m) {}
      SIG(const SIGModel::EpiPars* p) {
        nstates = SIGModel::nstates;
        pars = p;
        rho = p->rho;

        /* sampled and on ART (equivalent to recovery) **********************/

        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov",
                                                SIGModel::nstates,
                                                SIGModel::recoverChange,
                                                SIGModel::recovRateFun,
                                                SIGModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov",
                                              SIGModel::nstates,
                                              SIGModel::obsRecovChange,
                                              SIGModel::treeObsRecov,
                                              oneProb));
        simEvent.push_back(1);

        /* infections *******************************************************/

        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf",
                                                SIGModel::nstates,
                                                SIGModel::infChange,
                                                SIGModel::infRateFun,
                                                SIGModel::treeProbInf));

        obsTypes.push_back(new TransitionType("obsInf",
                                              SIGModel::nstates,
                                              SIGModel::obsInfChange,
                                              SIGModel::treeObsInf,
                                              oneProb));
        simEvent.push_back(1);

        /* sampled but not on ART *******************************************/

        typeMap[1] = 30;
        transTypes.push_back(new TransitionType("simInf",
                                                SIGModel::nstates,
                                                SIGModel::sampChange,
                                                SIGModel::sampRateFun,
                                                SIGModel::treeProbSample));

        obsTypes.push_back(new TransitionType("obsInf",
                                              SIGModel::nstates,
                                              SIGModel::obsSampleChange,
                                              SIGModel::treeObsSample,
                                              oneProb));
        simEvent.push_back(1);

      }

      virtual ~SIG() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      inline bool validState(const EpiState& es) const { return true; }
  };
}

#endif // __SIGMODEL_H__
