#include "Model.h"

namespace MCTraj {
  double Model::calculateTransRates(const EpiState& state, 
                                    vector<double>& transRates) const
  {
    size_t i = 0;
    transRates[i] = transTypes[i]->applyRate(state,pars);
    for (i = 1; i < transTypes.size(); ++i) {
      if (simEvent[i]) {
        transRates[i] = transRates[i-1] + transTypes[i]->applyRate(state,pars);
      } else {
        transRates[i] = transRates[i-1];
      }
    }
    return transRates[i-1];
  }

  size_t Model::chooseTransition(gsl_rng* rng, const vector<double>& transRates) const {
    double r = gsl_rng_uniform(rng)*transRates.back();
    int i = 0;
    int n = transRates.size();
    while (r > transRates[i] && i < n) ++i;
    return i;
  }

  bool Model::valid() const {
    if (transTypes.size() < nstates) return false;
    if (obsTypes.size() < nstates) return false;
    if (simEvent.size() < nstates) return false;
    return true;
  }
}

