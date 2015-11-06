#include "Parameters.h"

const parameter_t* Parameters::by_name(string name) const {
  map<string,int>::const_iterator i;
  i = name_map.find(name);
  if (i != name_map.end()) {
    return &pars[i->second];
  } else {
    return NULL;
  }
}

// ===========================================================================

size_t Parameters::nfree() const {
  size_t nf = 0;
  for (size_t i = 0; i < pars.size(); ++i) {
    if (! pars[i].lock) ++nf;
  }
  return nf;
}

// ===========================================================================

void Parameters::init_free_map() {
  free_map.clear();
  for (size_t i = 0; i < pars.size(); ++i) {
    if (! pars[i].lock) free_map.push_back(i);
  }
}

// ===========================================================================

double Parameters::limits(size_t j, double x, char trans) const{
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

// ===========================================================================

void Parameters::from_json(rapidjson::Document& jpars) {
  try {
    rapidjson::Value::MemberIterator m1;
    rapidjson::SizeType z; 

    m1 = jpars.FindMember("model_type");
    if (m1 != jpars.MemberEnd()) {
      if (! m1->value.IsString()) throw 101;
      model_type = m1->value.GetString();
    }

    m1 = jpars.FindMember("nroot");
    if (m1 != jpars.MemberEnd()) {
      if (! m1->value.IsInt()) throw 102;
      nroot = m1->value.GetInt();
    }

    m1 = jpars.FindMember("shifts");
    if (m1 != jpars.MemberEnd()) {
      if (! m1->value.IsArray()) throw 103;
      shifts.resize(m1->value.Size());
      for (rapidjson::SizeType i = 0; i < m1->value.Size(); ++i) {
        shifts[i] = m1->value[i].GetDouble();
      }
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

      m1 = a.FindMember("value"); 
      if (m1 == a.MemberEnd()) m1 = a.FindMember("init"); 
      if (m1 != a.MemberEnd()) {
        if (m1->value.IsDouble())  {
          pars[i].init.assign(1,m1->value.GetDouble());
        } else if (m1->value.IsArray()) {
          pars[i].init.assign(m1->value.Size(),0.0);
          for (size_t j = 0; j < m1->value.Size(); ++j) {
            pars[i].init[j] = m1->value[j].GetDouble();
          }
        } else {
          throw 2;
        }
      }

      m1 = a.FindMember("lock");
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsDouble()) throw 4;
        pars[i].lock = 1;
        pars[i].lo = m1->value.GetDouble();
        pars[i].hi = pars[i].lo;
        pars[i].init.assign(1,0.0);
        pars[i].init[0] = pars[i].hi;
      }

      m1 = a.FindMember("scale");
      if (m1 != a.MemberEnd()) {
        if (! m1->value.IsString()) throw 5;
        pars[i].scale = m1->value.GetString()[0];
      }
    }
  } catch (int e) {
    cerr << "Exception while reading pars: " << e << endl;
    return;
  }

  json_simpars(jpars);
}

// ===========================================================================

void Parameters::json_simpars(rapidjson::Document& jpars) {
  try {
    rapidjson::Value::MemberIterator it;

    it = jpars.FindMember("simulation");
    if (it != jpars.MemberEnd()) {
      if (! it->value.IsObject()) throw 201;
      rapidjson::Value& d = it->value;

      it = d.FindMember("maxTime");
      if (! it->value.IsDouble()) throw 202;
      sim_pars.maxTime = it->value.GetDouble();
    }
  } catch (int e) {
    cerr << "Exception while reading SimPars: " << e << endl;
    return;
  }
}

