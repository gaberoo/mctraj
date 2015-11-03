#ifndef __JSON_H__
#define __JSON_H__

#include <vector>

#include "Model.h"
#include "models/Models.h"

namespace MCTraj {
  string read_json(string filename);
  void get_json(rapidjson::Document& jpars, string json_input);
  void choose_model(Model*& mpt, EpiState*& es, PFPars& pf_pars, 
                    vector<Pars*>& vpars, const Parameters& p);
}

#endif
