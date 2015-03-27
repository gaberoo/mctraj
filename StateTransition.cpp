#include "StateTransition.h"

namespace MCTraj {
  ostream& operator<<(ostream& out, const StateTransition& trans) {
    out << trans.eventType << " " 
        << setw(12) << trans.time << " ";
    for (size_t i(0); i < trans.trans.size(); ++i) {
      out << setw(4) << trans.trans.at(i) << " ";
    }
    out << " ";
    for (size_t i(0); i < trans.branchTrans.size(); ++i) {
      out << trans.branchTrans.at(i).id << "_"
          << trans.branchTrans.at(i).new_color << " ";
    }
    return out;
  }

  string StateTransition::to_json() const {
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
    json_w.StartObject(); {
      json_w.String("time"); json_w.Double(time);
      json_w.String("prob"); json_w.Double(prob);
      json_w.String("trans"); json_w.StartArray(); {
        for (size_t i = 0; i < trans.size(); ++i) json_w.Int(trans[i]);
      } json_w.EndArray();
      json_w.String("type"); json_w.Int(eventType);
      json_w.String("absTime"); json_w.Double(absTime);
      json_w.String("branchChange"); json_w.StartArray();
      for (size_t i = 0; i < branchTrans.size(); ++i) {
        branchTrans[i].json(json_w);
      } json_w.EndArray();
    } json_w.EndObject();
    return buf.GetString();
  }
}

