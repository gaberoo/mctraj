#ifndef __MODEL_H__
#define __MODEL_H__

#include "MCTraj.h"
#include "EpiState.h"
#include "TransitionType.h"
#include "BranchState.h"

namespace MCTraj {
  inline double oneProb(const EpiState& es, const void* pars) { return 1.0; }

  class Model {
    public:
      Model() : nstates(0), pars(NULL), rho(0.0) {
        typeMap = new int[100];
        for (int i(0); i < 100; ++i) typeMap[i] = -1;
      }

      virtual ~Model() {
        delete[] typeMap;
      }

      size_t n() const { return nstates; }
      const void* p() const { return pars; }
      size_t ntrans() const { return transTypes.size(); }
      const TransitionType* getTType(size_t i) const { return transTypes[i]; }
      const TransitionType* getObsType(size_t i) const { return obsTypes[i]; }

      const TransitionType* getObsTypeMap(size_t i) const { 
        if (typeMap[i] >= 0) return obsTypes[typeMap[i]]; 
        else return NULL;
      }
      const TransitionType* getTTypeMap(size_t i) const { 
        if (typeMap[i] >= 0) return transTypes[typeMap[i]]; 
        else return NULL;
      }
      size_t mapType(size_t i) const { return typeMap[i]; }

      double calculateTransRates(const EpiState& state, vector<double>& transRates) const;

      inline const void* getPars() const { return pars; }
      inline void setPars(const void* p) { pars = p; }
      
      inline double getRho() const { return rho; }
      inline void setRho(double x) { rho = x; }

      size_t chooseTransition(gsl_rng* rng, const vector<double>& transRates) const;

      inline int& sim_event(size_t i) { return simEvent.at(i); }
      inline int sim_event(size_t i) const { return simEvent.at(i); }

      virtual double sample_rho(const EpiState& es, gsl_rng* rng, void* pars = NULL) const = 0;

      bool valid() const;

    protected:
      size_t nstates;
      const void* pars;
      double rho;
      int* typeMap;
      vector<const TransitionType*> transTypes;
      vector<const TransitionType*> obsTypes;
      vector<int> simEvent;
  };
}

#endif // __MODEL_H__
