#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <exception>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
using namespace std;

#include <rapidjson/document.h>

class parameter_t {
  public:
    parameter_t() 
      : name(""), init(1,0.0), lo(-INFINITY), hi(INFINITY), 
        lock(0), scale('n')
    {}
    virtual ~parameter_t() {}

    string name;
    vector<double> init;
    double lo;
    double hi;
    int lock;
    char scale;
};

typedef struct {
  double maxTime;
} SimPars;

class Parameters {
  public:
    Parameters() 
      : model_type("") 
    { 
      sim_pars.maxTime = 100.0; 
    }
    virtual ~Parameters() {}

    void from_json(rapidjson::Document&);

    const parameter_t* by_name(string name) const;
    size_t nfree() const;

    inline size_t nval(string name) const {
      const parameter_t* pt = by_name(name);
      return (pt != NULL) ? pt->init.size() : 0;
    }

    inline double value(string name, size_t pos = 0) const { 
      const parameter_t* pt = by_name(name);
      double x = 0.0;
      if (pt != NULL) {
        if (pt->init.size() == 0) return 0.0;
        else if (pos >= pt->init.size()) pos = pt->init.size()-1;
        switch (pt->scale) {
          case 'l':
          case 'L': x = exp(pt->init.at(pos)); break;
          default:  x = pt->init.at(pos);      break;
        }
      } else {
#ifdef DEBUG
        cerr << "Field '" << name << "' not found!" << endl;
#endif
      }
      return x;
    }

    void init_free_map();
    inline size_t free(size_t i) const { return free_map[i]; }
    inline double lower(size_t j) const { return pars[free_map[j]].lo; }
    inline double upper(size_t j) const { return pars[free_map[j]].hi; }
    double limits(size_t j, double x, char trans = 'n') const;

  protected:
    void json_simpars(rapidjson::Document& jpars);

  public:
    string model_type;
    int nroot;
    vector<double> shifts;
    vector<parameter_t> pars;
    map<string,int> name_map;
    vector<int> free_map;
    SimPars sim_pars;
};

// ===========================================================================

//class ParEx : public exception {
//  virtual const char* what() const throw() {
//    return "My exception happened";
//  }
//};

// ===========================================================================

#endif
