#include "Model.h"

namespace MCTraj {
  double Model::calculateTransRates(const EpiState& state, 
                                    vector<double>& transRates,
                                    vector<double>& trueRates) const
  {
    size_t i = 0;
    transRates[i] = transTypes[i]->applyRate(state,getPars(),trueRates[i]);
    for (i = 1; i < transTypes.size(); ++i) {
      if (simEvent[i]) {
        transRates[i] = transRates[i-1] + transTypes[i]->applyRate(state,getPars(),trueRates[i]);
      } else {
        trueRates[i] = 0.0;
        transRates[i] = transRates[i-1];
      }
    }
    return transRates[i-1];
  }

  /**************************************************************************/

  int Model::chooseTransition(rng::RngStream* rng, const vector<double>& transRates) const {
    return rng->pick(transRates.data(),transRates.size());
  }

  /**************************************************************************/

  double Model::delTransRate(vector<double>& transRates, size_t i) const {
    double rate = transRates[i] - ((i>0) ? transRates[i-1] : 0.0);
    for (size_t j = i; j < transRates.size(); ++j) transRates[j] -= rate;
    return rate;
  }

  /**************************************************************************/

  bool Model::valid() const {
    if (transTypes.size() < nstates) return false;
    if (obsTypes.size() < nstates) return false;
    if (simEvent.size() < nstates) return false;
    return true;
  }

  /**************************************************************************/

  string Pars::to_json() const {
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
    json(json_w);
    return buf.GetString();
  }
}

