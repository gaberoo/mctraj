#ifndef __SISMODEL_H__
#define __SISMODEL_H__

#include <rng/Rng.h>
#include "../Model.h"
#include "../EpiState.h"
#include "../Trajectory.h"

namespace MCTraj {
  namespace SISModel {
    class EpiPars : public Pars {
    public:
      EpiPars() {}
      EpiPars(const EpiPars& p) 
        : Pars(p), N(p.N), beta(p.beta), mu(p.mu), 
          psi(p.psi)
      {}
      virtual ~EpiPars() {}

      void json(rapidjson::Writer<rapidjson::StringBuffer>& w) const {
        w.StartObject(); {
          w.String("name");  w.String("SIS");
          w.String("N");     w.Double(N);
          w.String("beta");  w.Double(beta);
          w.String("mu");    w.Double(mu);
          w.String("psi");   w.Double(psi);
          w.String("rho");   w.Double(rho);
        } w.EndObject();
      }

      void from_parameters(const Parameters& p, size_t pos = 0) {
        N = p.value("N",pos);
        beta = p.value("beta",pos);
        mu = p.value("mu",pos);
        psi = p.value("psi",pos);
        rho = p.value("rho",pos);
      }

      double N;
      double beta;
      double mu;
      double psi;
    };

    /* S, I(all), I(tree) */
    const size_t nstates = 3;

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars, double& trueRate);
    double treeProbInf(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0 };

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars, double& trueRate);
    double treeProbRecov(const EpiState& es, const void* pars);
    const int recoverChange[] = { 1, -1, 0 };

    /************************************************************************/

    double treeObsInf(const EpiState& es, const void* pars, double& trueRate);
    const int obsInfChange[] = { -1, 1, 1 };

    /************************************************************************/

    double treeObsRecov(const EpiState& es, const void* pars, double& trueRate);
    const int obsRecovChange[] = { 1, -1, -1 };
  }

  /**************************************************************************/

  class SIS : public Model {
    public:
      SIS() {
        nstates = SISModel::nstates;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov", SISModel::nstates, SISModel::recoverChange, SISModel::recovRateFun, SISModel::treeProbRecov));
        obsTypes.push_back(new TransitionType("obsRecov", SISModel::nstates, SISModel::obsRecovChange, SISModel::treeObsRecov, oneProb));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf", SISModel::nstates, SISModel::infChange, SISModel::infRateFun, SISModel::treeProbInf));
        obsTypes.push_back(new TransitionType("obsInf", SISModel::nstates, SISModel::obsInfChange, SISModel::treeObsInf, oneProb));
        simEvent.push_back(1);
      }
      SIS(const SIS& m) : Model(m) {}

      virtual ~SIS() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      inline bool validState(const EpiState& es) const { return true; }
      void toTree(const Trajectory& traj, rng::RngStream* rng, vector<TreeNode>& tree) const;
  };
}

#endif // __SISMODEL_H__
