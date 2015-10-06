#ifndef __IMODEL_H__
#define __IMODEL_H__

#include <rng/Rng.h>
#include "../Model.h"
#include "../EpiState.h"
#include "SEIS.h"

namespace MCTraj {
  namespace IModel {
    typedef SEISModel::EpiPars EpiPars;

    /* I(all), I(tree) */
    const size_t nstates = 2;

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars, double& trueRate);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { 1, 0 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars, double& trueRate);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { -1, 0 };

    /************************************************************************/

    double treeObsInf(const EpiState& es, const void* pars, double& trueRate);
    const int obsInfChange[] = { 1, 1 };

    /************************************************************************/

    double treeObsRecov(const EpiState& es, const void* pars, double& trueRate);
    const int obsRecovChange[] = { -1, -1 };
  }

  /**************************************************************************/

  class I : public Model {
    public:
      I(const I& m) : Model(m) {}
      I(const IModel::EpiPars* p) {
        nstates = IModel::nstates;
        pars = p;
        rho = p->rho;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov",
                                                IModel::nstates,
                                                IModel::recoverChange,
                                                IModel::recovRateFun,
                                                IModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov",
                                              IModel::nstates,
                                              IModel::obsRecovChange,
                                              IModel::treeObsRecov,
                                              oneProb));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf",
                                                IModel::nstates,
                                                IModel::infChange,
                                                IModel::infRateFun,
                                                IModel::treeProbInf));

        obsTypes.push_back(new TransitionType("obsInf",
                                              IModel::nstates,
                                              IModel::obsInfChange,
                                              IModel::treeObsInf,
                                              oneProb));
        simEvent.push_back(1);
      }

      virtual ~I() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      inline bool validState(const EpiState& es) const { return true; }
      void toTree(const Trajectory& traj, rng::RngStream* rng, vector<TreeNode>& tree) const;
  };
}

#endif // __SISMODEL_H__
