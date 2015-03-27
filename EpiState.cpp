#include "EpiState.h"

namespace MCTraj {
  string EpiState::to_json() const {
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
    json(json_w);
    return buf.GetString();
  }
}

