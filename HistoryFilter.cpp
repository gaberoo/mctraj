#include "HistoryFilter.h"

namespace MCTraj {
  void MCTraj::HistoryFilter::setTree(const Tree* tree, int skip, rng::Rng* rng) {
      this->tree = tree; 
      // preallocate the space for the particles
      if (vflag > 1) cerr << "allocating..." << flush;
      pf.reserve(tree->size()+1);
      for (size_t i = 0; i <= tree->size(); ++i) {
        pf.push_back(vector<TrajParticle>(size()));
      }
      if (vflag > 1) cerr << "done." << endl;

      while (skip-- > 0) {
        addTreeEvent(model->p(),rng,1);

        inc();

        for (size_t i = 0; i < size(); ++i) {
          particle(i).copy(pf[curStep-1][i]);
          particle(i).setId(i);
          particle(i).setParent(i);
          particle(i).resetWeight();
          particle(i).setProb(1.0);
        }
      }

  }
}

