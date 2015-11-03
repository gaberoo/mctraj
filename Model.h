#ifndef __MODEL_H__
#define __MODEL_H__

#include <cstring>
#include <rng/Rng.h>

#include "MCTraj.h"
#include "EpiState.h"
#include "TransitionType.h"
#include "BranchState.h"
#include "Tree.h"
#include "TreeNode.h"
#include "Parameters.h"

namespace MCTraj {
  inline double oneProb(const EpiState& es, const void* pars) { return 1.0; }

  class Pars {
    public:
      Pars() : tree(NULL), rho(1.0) {}
      Pars(const Pars& p) : tree(p.tree), rho(p.rho) {}
      virtual ~Pars() {}

      virtual void json(rapidjson::Writer<rapidjson::StringBuffer>&) const {}
      string to_json() const;

      virtual void from_parameters(const Parameters& p, size_t pos = 0) = 0;

    public:
      const Tree* tree;
      double rho;
  };

  class Model {
    public:
      Model() 
        : nstates(0), pars_cnt(0), do_branches(false)
      {
        typeMap = new int[100];
        for (int i(0); i < 100; ++i) typeMap[i] = -1;
      }

      Model(const Model& m)
        : nstates(m.nstates), pars(m.pars), pars_cnt(m.pars_cnt), 
          transTypes(m.transTypes), obsTypes(m.obsTypes), 
          simEvent(m.simEvent), do_branches(m.do_branches)
      {
        typeMap = new int[100];
        memcpy(typeMap,m.typeMap,100*sizeof(int));
      }

      virtual ~Model() { delete[] typeMap; }

      size_t n() const { return nstates; }
      const void* p() const { 
        return (pars.size() > 0) ? pars.at(pars_cnt) : NULL; 
      }
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

      double calculateTransRates(const EpiState& state, 
                                 vector<double>& transRates,
                                 vector<double>& trueRates) const;
      double delTransRate(vector<double>& transRates, size_t i) const;

      inline const Pars* getPars() const { 
        return (pars.size() > 0) ? pars.at(pars_cnt) : NULL; 
      }
      inline void addPars(const Pars* p) { pars.push_back(p); }
      inline string pars_json(size_t pos = 0) const { 
        return (pars.size() > 0) ? pars.at(pos)->to_json() : ""; 
      }
      inline void incPars() {
        if ((int) pars.size() > pars_cnt+1) pars_cnt++;
      }
      
      inline double rho() const { 
        const Pars* p = getPars();
        return (p != NULL) ? p->rho : 1.0;
      }

      int chooseTransition(rng::RngStream* rng, const vector<double>& transRates) const;

      inline int& sim_event(size_t i) { return simEvent.at(i); }
      inline int sim_event(size_t i) const { return simEvent.at(i); }

      bool valid() const;

      virtual double sample_rho(const EpiState& es, rng::RngStream* rng, void* pars = NULL) const = 0;
      virtual bool validState(const EpiState& es) const = 0;

      inline size_t nTransRates() const { return ntrans(); }

      virtual void toTree(const Trajectory& traj, rng::RngStream* rng, vector<TreeNode>& tree) const = 0;

    protected:
      size_t nstates;
      vector<const Pars*> pars;
      int pars_cnt;
      int* typeMap;
      vector<const TransitionType*> transTypes;
      vector<const TransitionType*> obsTypes;
      vector<int> simEvent;

    public:
      bool do_branches;
  };
}

#endif // __MODEL_H__
