#include "MCTraj.h"

void MCTraj::pf_pars_read_json(PFPars* p, rapidjson::Value& json) {
  // reading priors from JSON
  rapidjson::Value::MemberIterator m1;

  // find variable ranges
  rapidjson::Value::MemberIterator it = json.FindMember("pf");
  if (it == json.MemberEnd()) throw "No pf section in JSON.";

  rapidjson::Value& d = it->value;

  p->reps = 1;
  p->full_tree = 0;

  m1 = d.FindMember("model_type");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsString()) throw "model_type";
    p->model_type = m1->value.GetString()[0];
  }

  m1 = d.FindMember("num_particles");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "num_particles";
    p->num_particles = m1->value.GetInt();
  }

  m1 = d.FindMember("print_particles");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "print_particles";
    p->print_particles = m1->value.GetInt();
  }

  m1 = d.FindMember("skip");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "skip";
    p->skip = m1->value.GetInt();
  }

  m1 = d.FindMember("save_traj");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "save_traj";
    p->print_traj = m1->value.GetInt();
  }

  m1 = d.FindMember("seed");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "seed";
    p->seed = m1->value.GetInt();
  }

  m1 = d.FindMember("filter_time");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsDouble()) throw "filter_time";
    p->filter_time = m1->value.GetDouble();
  }

  m1 = d.FindMember("step_size");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsDouble()) throw "step_size";
    p->step_size = m1->value.GetDouble();
  }

  m1 = d.FindMember("adj_zero");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "adj_zero";
    p->adj_zero = m1->value.GetInt();
  }

  m1 = d.FindMember("vflag");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "vflag";
    p->vflag = m1->value.GetInt();
  }

  m1 = d.FindMember("history");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "history";
    p->history = m1->value.GetInt();
  }
}

