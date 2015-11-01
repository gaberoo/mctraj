#ifndef __PSO_PARS_H__
#define __PSO_PARS_H__

#include <iostream>
#include <rapidjson/document.h>

class pso_pars_t {
public:
  pso_pars_t()
    : max_evals(100), num_particles(20), init_type(1),
      vflag(0), slowdown(0), out(NULL), hist(NULL)
  {}

  virtual ~pso_pars_t() {
    if (out != NULL) delete out;
    if (hist != NULL) delete hist;
  }

  void from_json(rapidjson::Value& d);

  int max_evals;
  int num_particles;
  int init_type;
  int vflag;
  int slowdown;
  std::ostream* out;
  std::ostream* hist;
};

inline void pso_pars_t::from_json(rapidjson::Value& d) 
{
  rapidjson::Value::MemberIterator it;

  it = d.FindMember("max_evals");
  if (it != d.MemberEnd()) max_evals = it->value.GetInt();

  it = d.FindMember("num_particles");
  if (it != d.MemberEnd()) num_particles = it->value.GetInt();

  it = d.FindMember("init_type");
  if (it != d.MemberEnd()) init_type = it->value.GetInt();

  it = d.FindMember("vflag");
  if (it != d.MemberEnd()) vflag = it->value.GetInt();

  it = d.FindMember("slowdown");
  if (it != d.MemberEnd()) slowdown = it->value.GetInt();

  it = d.FindMember("out_fn");
  if (it != d.MemberEnd()) out = new ofstream(it->value.GetString());

  it = d.FindMember("hist_fn");
  if (it != d.MemberEnd()) hist = new ofstream(it->value.GetString());
}


#endif

