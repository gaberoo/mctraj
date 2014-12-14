#include "TransitionType.h"

namespace MCTraj {
  ostream& operator<<(ostream& out, const TransitionType& tt) {
    for (size_t i(0); i < tt.change.size(); ++i) {
      out << tt.change[i] << " ";
    }
    return out;
  }

}
