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
        : Pars(p), N(p.N), beta(p.beta), mu(p.mu), psi(p.psi)
      {}
      virtual ~EpiPars() {}

      inline void json(rapidjson::Writer<rapidjson::StringBuffer>& w) const {
        w.StartObject(); {
          w.String("name");  w.String("SIS");
          w.String("N");     w.Double(N);
          w.String("beta");  w.Double(beta);
          w.String("mu");    w.Double(mu);
          w.String("psi");   w.Double(psi);
          w.String("rho");   w.Double(rho);
        } w.EndObject();
      }

      inline void from_parameters(const Parameters& p, size_t pos = 0) {
        N = p.value("N",pos);
        beta = p.value("beta",pos);
        mu = p.value("mu",pos);
        psi = p.value("psi",pos);
        rho = p.value("rho",pos);
      }

      inline void from_state(const double* state, const char* scales) {
        N     = (scales[0] == 'l') ? exp(state[0]) : state[0];
        beta  = (scales[1] == 'l') ? exp(state[1]) : state[1];
        mu    = (scales[2] == 'l') ? exp(state[2]) : state[2];
        psi   = (scales[3] == 'l') ? exp(state[3]) : state[3];
        rho   = (scales[4] == 'l') ? exp(state[4]) : state[4];
      }

      inline void to_state(double* state, const char* scales) {
        cerr << scales << endl;
        state[0] = (scales[0] == 'l') ? log(N)    : N;
        state[1] = (scales[1] == 'l') ? log(beta) : beta;
        state[2] = (scales[2] == 'l') ? log(mu)   : mu;
        state[3] = (scales[3] == 'l') ? log(psi)  : psi;
        state[4] = (scales[4] == 'l') ? log(rho)  : rho;
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

      inline void addPars(Pars* p) {
        try {
          SISModel::EpiPars* pp = dynamic_cast<SISModel::EpiPars*>(p);
          pars.push_back(new SISModel::EpiPars(*pp));
        } catch (exception& e) {
          cerr << "Error casting pointer!" << endl;
          abort();
        }
      }

      inline EpiState initState(const Parameters& p) const {
        SISModel::EpiPars epi = *dynamic_cast<SISModel::EpiPars*>(pars.front());
        EpiState init(SISModel::nstates);
        init[0] = (int) epi.N - p.nroot;
        init[1] = p.nroot;
        init[2] = p.nroot;
        return init;
      }
  };
}

#endif // __SISMODEL_H__
