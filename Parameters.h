#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <exception>
#include <rapidjson/document.h>
#include <vector>
#include <map>
using namespace std;

typedef struct {
  string name;
  double init;
  double lo;
  double hi;
  int lock;
} parameter_t;

class Parameters {
public:
  Parameters() : model_type("") {}
  virtual ~Parameters() {}

  void from_json(rapidjson::Document&);

  const parameter_t* by_name(string name) const;
  size_t nfree() const;

  void init_free_map();
  inline size_t free(size_t i) const { return free_map[i]; }
  inline double lower(size_t j) const { return pars[free_map[j]].lo; }
  inline double upper(size_t j) const { return pars[free_map[j]].hi; }
  double limits(size_t j, double x, char trans = 'n') const;

  string model_type;
  vector<parameter_t> pars;
  map<string,int> name_map;
  vector<int> free_map;
};

//class ParEx : public exception {
//  virtual const char* what() const throw() {
//    return "My exception happened";
//  }
//};

inline const parameter_t* Parameters::by_name(string name) const {
  map<string,int>::const_iterator i;
  i = name_map.find(name);
  if (i != name_map.end()) {
    return &pars[i->second];
  } else {
    return NULL;
  }
}

inline size_t Parameters::nfree() const {
  size_t nf = 0;
  for (size_t i = 0; i < pars.size(); ++i) {
    if (! pars[i].lock) ++nf;
  }
  return nf;
}

inline void Parameters::init_free_map() {
  free_map.clear();
  for (size_t i = 0; i < pars.size(); ++i) {
    if (! pars[i].lock) free_map.push_back(i);
  }
}

inline double Parameters::limits(size_t j, double x, char trans) const{
  double y = x;
  switch (trans) {
    case 'l':
    case 'L':
      if (x < log(lower(j))) y = log(lower(j));
      else if (x > log(upper(j))) y = log(upper(j));
      break;
    case 'n':
    case 'N':
    default:
      if (x < lower(j)) y = lower(j);
      else if (x > upper(j)) y = upper(j);
      break;
  }
  return y;
}

inline void Parameters::from_json(rapidjson::Document& jpars) {
  try {
    rapidjson::Value::MemberIterator m1;
    rapidjson::SizeType z; 

    m1 = jpars.FindMember("model_type");
    if (m1 != jpars.MemberEnd()) {
      if (! m1->value.IsString()) throw 101;
      model_type = m1->value.GetString();
    }

    rapidjson::Value::MemberIterator _d = jpars.FindMember("pars");
    if (_d == jpars.MemberEnd()) throw 1;

    rapidjson::Value& d = _d->value;

    pars.resize(d.Size());

    for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
      rapidjson::Value& a = d[i];

      m1 = a.FindMember("name"); 
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsString()) throw 5;
        pars[i].name = m1->value.GetString();
        name_map.insert(make_pair(pars[i].name,i));
      }

      m1 = a.FindMember("limits"); 
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsArray()) throw 2;
        if (m1->value.Size() != 2) throw 3;
        z = 0; pars[i].lo = m1->value[z].GetDouble();
        z = 1; pars[i].hi = m1->value[z].GetDouble();
      }

      m1 = a.FindMember("init"); 
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsDouble()) throw 2;
        pars[i].init = m1->value.GetDouble();
      }

      m1 = a.FindMember("lock");
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsDouble()) throw 4;
        pars[i].lock = 1;
        pars[i].init = pars[i].lo = pars[i].hi = m1->value.GetDouble();
      }
    }
  } catch (int e) {
    cerr << "Exception while reading pars: " << e << endl;
    return;
  }
}


#endif
