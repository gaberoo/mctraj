#include "StaticFilter.h"

namespace MCTraj {
  StaticFilter::~StaticFilter() {
//    size_t i;
//#pragma omp parallel for
//    for (i = 0; i < pf1.size(); ++i) { pf1[i].~TrajParticle(); }
//#pragma omp parallel for
//    for (i = 0; i < pf2.size(); ++i) { pf2[i].~TrajParticle(); }
  }

  // =========================================================================

  void StaticFilter::meanState(double x[], size_t gen) const {
    size_t i, j;
    size_t n = model->n();
    EpiState y;
    for (i = 0; i < n; ++i) x[i] = 0.0;
    for (j = 0; j < size(); ++j) {
      y = particle(j).getState();
      for (i = 0; i < n; ++i) x[i] += y[i];
    }
    for (i = 0; i < n; ++i) x[i] /= 1.0*size();
  }

  // =========================================================================

  void StaticFilter::varState(double var[], const double mean[], size_t gen) const {
    size_t i, j;
    size_t n = model->n();
    EpiState y;
    double z;
    for (i = 0; i < n; ++i) var[i] = 0.0;
    for (j = 0; j < size(); ++j) {
      y = particle(j).getState();
      for (i = 0; i < n; ++i) {
        z = y[i]-mean[i];
        var[i] += z*z;
      }
    }
    for (i = 0; i < n; ++i) var[i] /= 1.0*(size()-1);
  }

  // =========================================================================

  void StaticFilter::copyFromPrev(size_t i, size_t j) {
    (*A)[i].setState((*B)[j].getState());
  }

  // =========================================================================

  Trajectory StaticFilter::singleTraj(int p) const {
    return particle(p);
  }
}
