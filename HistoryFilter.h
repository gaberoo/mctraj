#ifndef __HISTORYFILTER_H__
#define __HISTORYFILTER_H__

#include "TrajParticleFilter.h"

namespace MCTraj {
  class HistoryFilter : public TrajParticleFilter {
    public:
      HistoryFilter(const Model* m) : TrajParticleFilter(m) {}
      HistoryFilter(const HistoryFilter& x) : TrajParticleFilter(x) {}
      virtual ~HistoryFilter() {}

      void setTree(const Tree* tree, int skip = 0, rng::Rng* rng = NULL);
  };
}

#endif
