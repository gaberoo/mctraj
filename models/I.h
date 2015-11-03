#ifndef __IMODEL_H__
#define __IMODEL_H__

#include <rng/RngStream.h>
#include <gsl/gsl_math.h>

#include "../Model.h"
#include "../EpiState.h"
#include "../StateTransition.h"
#include "../ascii.h"

namespace MCTraj {
  namespace IModel {
    class EpiPars : public Pars {
    public:
      EpiPars() {}
      EpiPars(const EpiPars& e) 
        : Pars(e), beta(e.beta), mu(e.mu), psi(e.psi)
      {}
      virtual ~EpiPars() {}

      void json(rapidjson::Writer<rapidjson::StringBuffer>& w) const {
        w.StartObject(); {
          w.String("name");  w.String("I");
          w.String("beta");  w.Double(beta);
          w.String("mu");    w.Double(mu);
          w.String("psi");   w.Double(psi);
          w.String("rho");   w.Double(rho);
        } w.EndObject();
      }

      void from_parameters(const Parameters& p, size_t pos = 0) {
        beta = p.value("beta",pos);
        mu = p.value("mu",pos);
        psi = p.value("psi",pos);
        rho = p.value("rho",pos);
      }

      double beta;   // per-contact infection rate
      double mu;     // recovery rate
      double psi;    // sampling rate
    };

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
      I() {
        nstates = IModel::nstates;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov", IModel::nstates, IModel::recoverChange, IModel::recovRateFun, IModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov", IModel::nstates, IModel::obsRecovChange, IModel::treeObsRecov, oneProb));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf", IModel::nstates, IModel::infChange, IModel::infRateFun, IModel::treeProbInf));
        obsTypes.push_back(new TransitionType("obsInf", IModel::nstates, IModel::obsInfChange, IModel::treeObsInf, oneProb));
        simEvent.push_back(1);
      }
      I(const I& m) : Model(m) {}
      virtual ~I() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      inline bool validState(const EpiState& es) const { return true; }
      void toTree(const Trajectory& traj, rng::RngStream* rng, vector<TreeNode>& tree) const;
  };
}

#endif // __SISMODEL_H__
