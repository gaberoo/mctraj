#include "MCTraj.h"

void MCTraj::PFPars::from_json(rapidjson::Value& json) {
  // reading priors from JSON
  rapidjson::Value::MemberIterator m1;

  // find variable ranges
  rapidjson::Value::MemberIterator it = json.FindMember("pf");
  if (it == json.MemberEnd()) throw "No pf section in JSON.";

  rapidjson::Value& d = it->value;

  full_tree = 0;

  m1 = d.FindMember("model_type");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsString()) throw "model_type";
    model_type = m1->value.GetString()[0];
  }

  m1 = d.FindMember("num_particles");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "num_particles";
    num_particles = m1->value.GetInt();
  }

  m1 = d.FindMember("print_particles");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "print_particles";
    print_particles = m1->value.GetInt();
  }

  m1 = d.FindMember("skip");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "skip";
    skip = m1->value.GetInt();
  }

  m1 = d.FindMember("save_traj");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "save_traj";
    print_traj = m1->value.GetInt();
  }

  m1 = d.FindMember("seed");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "seed";
    seed = m1->value.GetInt();
  }

  m1 = d.FindMember("filter_time");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsDouble()) throw "filter_time";
    filter_time = m1->value.GetDouble();
  }

  m1 = d.FindMember("step_size");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsDouble()) throw "step_size";
    step_size = m1->value.GetDouble();
  }

  m1 = d.FindMember("adj_zero");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "adj_zero";
    adj_zero = m1->value.GetInt();
  }

  m1 = d.FindMember("vflag");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "vflag";
    vflag = m1->value.GetInt();
  }

  m1 = d.FindMember("history");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "history";
    history = m1->value.GetInt();
  }

  m1 = d.FindMember("reps");
  if (m1 != d.MemberEnd()) {
    if (! m1->value.IsInt()) throw "reps";
    reps = m1->value.GetInt();
  }

}


