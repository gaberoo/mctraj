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
}

