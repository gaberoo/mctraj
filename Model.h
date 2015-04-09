#ifndef __MODEL_H__
#define __MODEL_H__

#include <cstring>
#include <rng/Rng.h>

#include "MCTraj.h"
#include "EpiState.h"
#include "TransitionType.h"
#include "BranchState.h"

namespace MCTraj {
  inline double oneProb(const EpiState& es, const void* pars) { return 1.0; }

  class Pars {
    public:
      Pars() {}
      virtual ~Pars() {}

      virtual void json(rapidjson::Writer<rapidjson::StringBuffer>&) const {}
      string to_json() const;
  };

  class Model {
    public:
      Model() : nstates(0), pars(NULL), rho(0.0) {
        typeMap = new int[100];
        for (int i(0); i < 100; ++i) typeMap[i] = -1;
      }

      Model(const Model& m)
        : nstates(m.nstates), pars(m.pars), rho(m.rho),
          transTypes(m.transTypes), obsTypes(m.obsTypes), 
          simEvent(m.simEvent)
      {
        typeMap = new int[100];
        memcpy(typeMap,m.typeMap,100*sizeof(int));
      }

      virtual ~Model() { delete[] typeMap; }

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
      double delTransRate(vector<double>& transRates, size_t i) const;

      inline const Pars* getPars() const { return pars; }
      inline void setPars(const Pars* p) { pars = p; }
      inline string pars_json() const { return pars->to_json(); }
      
      inline double getRho() const { return rho; }
      inline void setRho(double x) { rho = x; }

      size_t chooseTransition(rng::RngStream* rng, const vector<double>& transRates) const;

      inline int& sim_event(size_t i) { return simEvent.at(i); }
      inline int sim_event(size_t i) const { return simEvent.at(i); }

      bool valid() const;

      virtual double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const = 0;
      virtual bool validState(const EpiState& es) const = 0;

    protected:
      size_t nstates;
      // const void* pars;
      const Pars* pars;
      double rho;
      int* typeMap;
      vector<const TransitionType*> transTypes;
      vector<const TransitionType*> obsTypes;
      vector<int> simEvent;
  };
}

#endif // __MODEL_H__
