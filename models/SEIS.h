#ifndef __SEISMODEL_H__
#define __SEISMODEL_H__

#include <rng/RngStream.h>

#include <gsl/gsl_math.h>
#include "../Model.h"
#include "../EpiState.h"
#include "../StateTransition.h"

#include "../ascii.h"

namespace MCTraj {
  namespace SEISModel {
    class EpiPars : public Pars {
    public:
      EpiPars() : alpha(10.0) {}
      EpiPars(const EpiPars& e) 
        : Pars(e), N(e.N), beta(e.beta), mu(e.mu), 
          psi(e.psi), gamma(e.gamma), alpha(e.alpha)
      {}
      virtual ~EpiPars() {}

      void json(rapidjson::Writer<rapidjson::StringBuffer>& w) const {
        w.StartObject(); {
          w.String("name");  w.String("SEIS");
          w.String("N");     w.Double(N);
          w.String("beta");  w.Double(beta);
          w.String("mu");    w.Double(mu);
          w.String("psi");   w.Double(psi);
          w.String("rho");   w.Double(rho);
          w.String("gamma"); w.Double(gamma);
        } w.EndObject();
      }

      void from_parameters(const Parameters& p, size_t pos = 0) {
        N = p.value("N",pos);
        beta = p.value("beta",pos);
        mu = p.value("mu",pos);
        psi = p.value("psi",pos);
        rho = p.value("rho",pos);
        gamma = p.value("gamma",pos);
        alpha = p.value("alpha",pos);
      }

      inline void from_state(const double* state, const char* scales) {
        N     = (scales[0] == 'l') ? exp(state[0]) : state[0];
        beta  = (scales[1] == 'l') ? exp(state[1]) : state[1];
        mu    = (scales[2] == 'l') ? exp(state[2]) : state[2];
        psi   = (scales[3] == 'l') ? exp(state[3]) : state[3];
        rho   = (scales[4] == 'l') ? exp(state[4]) : state[4];
        gamma = (scales[5] == 'l') ? exp(state[5]) : state[5];
      }

      inline void to_state(double* state, const char* scales) {
        state[0] = (scales[0] == 'l') ? log(N)     : N;
        state[1] = (scales[1] == 'l') ? log(beta)  : beta;
        state[2] = (scales[2] == 'l') ? log(mu)    : mu;
        state[3] = (scales[3] == 'l') ? log(psi)   : psi;
        state[4] = (scales[4] == 'l') ? log(rho)   : rho;
        state[5] = (scales[5] == 'l') ? log(gamma) : gamma;
      }

      double N;      // total population size
      double beta;   // per-contact infection rate
      double mu;     // recovery rate
      double psi;    // sampling rate
      double gamma;  // transition rate

      double alpha;  // damping coefficient of the potential
    };

    const size_t nstates = 3; /* S, E, I, kE, kI */

    /************************************************************************/

    double infRateFun(const EpiState& es, const void* pars, double& trueRate);
    double infTreeProb(const EpiState& es, const void* pars);
    const int infChange[] = { -1, 1, 0 };
    int infValidate(const EpiState& es, const void* pars);

    /************************************************************************/

    double transRateFun(const EpiState& es, const void* pars, double& trueRate);
    double transTreeProb(const EpiState& es, const void* pars);
    const int transChange[] = { 0, -1, 1 };
    int transValidate(const EpiState& es, const void* pars);

    /************************************************************************/

    double recovRateFun(const EpiState& es, const void* pars, double& trueRate);
    double recovTreeProb(const EpiState& es, const void* pars);
    const int recoverChange[] = { 1, 0, -1 };
    int recovValidate(const EpiState& es, const void* pars);

    /************************************************************************/

    double infTreeObs(const EpiState& es, const void* pars, double& trueRate);
    const int infObsChange[] = { -1, 1, 0 };

    /************************************************************************/

    double transTreeObs(const EpiState& es, const void* pars, double& trueRate);
    const int obsTransChange[] = { 0, -1, 1 };

    /************************************************************************/

    double recovTreeObs(const EpiState& es, const void* pars, double& trueRate);
    const int recovChangeObs[] = { 1, 0, -1 };

    /************************************************************************/

    int infBranch(const EpiState& es, rng::RngStream* rng, 
                  StateTransition& st, const void* pars);
    int recovBranch(const EpiState& es, rng::RngStream* rng, 
                    StateTransition& st, const void* pars);
    int transBranch(const EpiState& es, rng::RngStream* rng, 
                    StateTransition& st, const void* pars);

    int infBranchObs(const EpiState& es, rng::RngStream* rng, 
                     StateTransition& st, const void* pars);
    int recovBranchObs(const EpiState& es, rng::RngStream* rng, 
                       StateTransition& st, const void* pars);
    int transBranchObs(const EpiState& es, rng::RngStream* rng, 
                       StateTransition& st, const void* pars);

    /************************************************************************/

    double calcBranchPotentials
      (const EpiState& es, const void* pars, int color, double* rates = NULL);

    const int branchTransChange[] = { -1, 1 };
    const int branchTransNew[] = { 0, 1 };
    const int branchTransOld[] = { 1, 0 };
  }

  /**************************************************************************/

  class SEIS : public Model {
    public:
      SEIS() {
        nstates = SEISModel::nstates;
        do_branches = true;

        /* recovery events */
        typeMap[0] = 0;
        transTypes.push_back(new TransitionType("simRecov", SEISModel::nstates, SEISModel::recoverChange, SEISModel::recovRateFun, SEISModel::recovTreeProb, SEISModel::recovBranch));
        obsTypes.push_back(new TransitionType("obsRecov", SEISModel::nstates, SEISModel::recovChangeObs, SEISModel::recovTreeObs, oneProb, SEISModel::recovBranchObs));
        simEvent.push_back(1);

        /* infection events */
        typeMap[1] = 1;
        transTypes.push_back(new TransitionType("simInf", SEISModel::nstates, SEISModel::infChange, SEISModel::infRateFun, SEISModel::infTreeProb, SEISModel::infBranch));
        obsTypes.push_back(new TransitionType("obsInf", SEISModel::nstates, SEISModel::infObsChange, SEISModel::infTreeObs, oneProb, SEISModel::infBranchObs));
        simEvent.push_back(1);

        /* transition events */
        typeMap[70] = 2;
        transTypes.push_back(new TransitionType("simTrans", SEISModel::nstates, SEISModel::transChange, SEISModel::transRateFun, SEISModel::transTreeProb, SEISModel::transBranch));
        obsTypes.push_back(new TransitionType("obsTrans", SEISModel::nstates, SEISModel::obsTransChange, SEISModel::transTreeObs, oneProb, SEISModel::transBranchObs));
        simEvent.push_back(1);

      }
      SEIS(const SEIS& m) : Model(m) { do_branches = true; }
      virtual ~SEIS() {}

      double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const;
      bool validState(const EpiState& es) const;
      void toTree(const Trajectory& traj, rng::RngStream* rng, vector<TreeNode>& tree) const;

      inline void addPars(Pars* p) {
        try {
          SEISModel::EpiPars* pp = dynamic_cast<SEISModel::EpiPars*>(p);
          pars.push_back(new SEISModel::EpiPars(*pp));
        } catch (exception& e) {
          cerr << "Error casting pointer!" << endl;
          abort();
        }
      }

      inline EpiState initState(const Parameters& p) const {
        SEISModel::EpiPars epi = *dynamic_cast<SEISModel::EpiPars*>(pars.front());
        EpiState init(SEISModel::nstates);
        init[0] = (int) epi.N - 1;
        init[1] = 1;
        init[2] = 0;
        return init;
      }
  };
}

#endif // __SEISMODEL_H__
