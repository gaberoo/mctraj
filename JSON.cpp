#include "JSON.h"

namespace MCTraj {
  // =========================================================================
  // read JSON file

  string read_json(string json_fn) {
    ifstream in;
    string json_input;
    if (json_fn != "") {
      in.open(json_fn.c_str());
      in.seekg(0,ios::end);
      json_input.reserve(in.tellg());
      in.seekg(0,ios::beg);
      json_input.assign(istreambuf_iterator<char>(in), 
                        istreambuf_iterator<char>());
    }
    return json_input;
  }

  // =========================================================================
  // Extract parameters from JSON

  void get_json(rapidjson::Document& jpars, string json_input) {
    try {
      jpars.Parse<0>(json_input.c_str());
      if (! jpars.IsObject()) throw 10;
    } catch (int e) {
      cerr << "Coudln't create JSON document:" << endl;
      cerr << json_input << endl;
    }
  }

  // =========================================================================
  // Extract parameters from JSON

  void choose_model(Model*& mpt, 
                    EpiState*& es, 
                    PFPars& pf_pars,
                    vector<Pars*>& vpars,
                    const Parameters& p)
  {
    size_t npars = p.shifts.size()+1;
    size_t ip = 0;

    SISModel::EpiPars  *sis_pars = NULL;
    SEISModel::EpiPars *seis_pars = NULL;

    switch (pf_pars.model_type) {
      case 'I':
      case 0:
        if (pf_pars.vflag > 1) cerr << "I model." << endl;
        mpt = new I;
        for (ip = 0; ip < npars; ++ip) {
          vpars.push_back(new IModel::EpiPars);
          vpars.back()->from_parameters(p,ip);
          mpt->addPars(vpars.back());
        }
        es = new EpiState(IModel::nstates);
        (*es)[0] = 0;
        break;

      case 'S':
      case 1:
      default:
        if (pf_pars.vflag > 1) cerr << "SIS model." << endl;
        mpt = new SIS;
        for (ip = 0; ip < npars; ++ip) {
          vpars.push_back(new SISModel::EpiPars);
          vpars.back()->from_parameters(p,ip);
          mpt->addPars(vpars.back());
        }
        es = new EpiState(SISModel::nstates);
        sis_pars = (SISModel::EpiPars*) vpars.front();
        (*es)[0] = (int) sis_pars[0].N - p.nroot;
        (*es)[1] = p.nroot;
        (*es)[2] = p.nroot;
        break;

      case 'R':
      case 2:
        if (pf_pars.vflag > 1) cerr << "SIR model." << endl;
        mpt = new SIR;
        for (ip = 0; ip < npars; ++ip) {
          vpars.push_back(new SISModel::EpiPars);
          vpars.back()->from_parameters(p,ip);
          mpt->addPars(vpars.back());
        }
        es = new EpiState(SISModel::nstates);
        sis_pars = (SISModel::EpiPars*) vpars.front();
        (*es)[0] = (int) sis_pars[0].N - p.nroot;
        (*es)[1] = p.nroot;
        (*es)[2] = p.nroot;
        break;

      case 'E':
      case 3:
        if (pf_pars.vflag > 1) cerr << "SEIS model." << endl;
        mpt = new SEIS;
        for (ip = 0; ip < npars; ++ip) {
          vpars.push_back(new SEISModel::EpiPars);
          vpars.back()->from_parameters(p,ip);
          mpt->addPars(vpars.back());
        }
        es = new EpiState(SEISModel::nstates);
        seis_pars = (SEISModel::EpiPars*) vpars.front();
        (*es)[0] = ((int) seis_pars[0].N)-1;
        (*es)[1] = 1;
        (*es)[2] = 0;
        if (pf_pars.full_tree) {
          mpt->sim_event(0) = 0;
          mpt->sim_event(1) = 0;
        }
        break;
    }
  }
}


