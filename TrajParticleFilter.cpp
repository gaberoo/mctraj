#include "TrajParticleFilter.h"

namespace MCTraj {
  void TrajParticleFilter::setTree(const Tree* tree, int skip, rng::Rng* rng) 
  { 
    this->tree = tree; 
  }

  // =========================================================================

  int TrajParticleFilter::filter(rng::Rng* rng, char filter_type) 
  {
    size_t n(size());

    // calculate total weight
    vector<double> w(n,0.0);
    double totalW = w_vec(w);
    if (vflag > 1) {
      double s2 = gsl_stats_variance(w.data(),1,w.size());
      cerr << "  FIL (" << filter_type << ") :: Sum = " << totalW 
           << ", var(w) = " << sqrt(s2)/(totalW/n) << endl;
      size_t k = gsl_stats_max_index(w.data(),1,n);
      cerr << "  Best particle = " << particle(k).getState() << endl;
    }

    if (totalW <= 0.0) { logw = -INFINITY; return -1; }

    // allocate space for new generation of particles
    if (vflag > 1) cerr << "  Increasing array size..." << flush;
    inc();
    if (vflag > 1) cerr << "done." << endl;

    if (filter_type == 'c') {
      // copy filter
      size_t i;
#pragma omp parallel for default(shared) private(i)
      for (i = 0; i < n; ++i) {
        particle(i).copy(prev(i));
        particle(i).setId(i);
        particle(i).setParent(i);
      }
    } else {
      // bootstrap filter
      logw += log(totalW/n);

      size_t i;

      if (vflag > 1) cerr << "  Preparing weight vector..." << flush;
      vector<double> x(n+1,0.0);
      rng::make_discrete(n,w.data(),x.data());
      if (vflag > 1) cerr << "done." << endl;

      size_t id;
      int p;

#pragma omp parallel for default(shared) private(i,id,p)
      for (i = 0; i < n; ++i) {
        try {
          // sampling with replacement
          id = omp_get_thread_num();
          p = (*rng)[id]->discrete_x(n,x.data());
//          if (w[p] <= 0.0) {
//            cerr << " !!! Filtered a particle with weight zero!" 
//                 << " p = " << p << "; ";
//            size_t start = r/x.back()*n;
//            cerr << " start = " << start << " "; 
//            if (p > 0) cerr << "(p-1) = " << x[p-1] << " ";
//            else cerr << "--- ";
//            cerr << "p = " << x[p] << " ";
//            if (p < n-1) cerr << "(p+1) = " << x[p+1] << " ";
//            else cerr << "--- ";
//            cerr << x.back() << " " << r << endl;
//          }
          particle(i).copy(prev(p));
          particle(i).setId(i);
          particle(i).setParent(p);
          particle(i).resetWeight();
          particle(i).resetProb();
        } catch (std::exception& e) {
          cerr << "Exception caught: " << e.what() << endl;
          abort();
        }
      }
    }

    if (vflag > 1) cerr << "Done filtering." << endl;

    return 1;
  }

  // ========================================================================

  double TrajParticleFilter::w_vec(vector<double>& w) const {
    size_t n(size());
    double totalW(0.0);
    w.resize(n);
    size_t i;
    for (i = 0; i < n; ++i) {
      w[i] = particle(i).getWeight();
      totalW += w[i];
    }
    return totalW;
  }

  // ========================================================================

  double TrajParticleFilter::p_vec(vector<double>& p) const {
    size_t n(size());
    double totalP(0.0);
    p.resize(n);
    size_t i;
    for (i = 0; i < n; ++i) {
      p[i] = particle(i).getProb();
      totalP += p[i];
    }
    return totalP;
  }

  // ========================================================================

  TrajParticle TrajParticleFilter::sample(rng::RngStream* rng) const {
    size_t n(size());
    vector<double> w(n);
    w_vec(w);
    size_t j = rng->discrete(n,w.data());
    return particle(j);
  }

  // ========================================================================

  int TrajParticleFilter::addTreeEvent(const void* pars, rng::Rng* rng, int noProb) {
    size_t nextStep = curStep;

    double eventTime = tree->times[nextStep];
    int eventType = tree->ttypes[nextStep];
    int modelType = model->mapType(eventType);

    if (vflag > 1) {
      cerr << "Add tree event :: " << setw(12) << eventTime << " " 
           << "(" << eventType << "->" << modelType << ") :: ";
      for (size_t i(0); i < tree->ids[nextStep].size(); ++i) {
        cerr << tree->ids[nextStep][i] << " ";
      }
      cerr << endl;
    }
    if (modelType < 0) return tree->ttypes[nextStep];

    size_t j;
    int id;
    double dw;

#pragma omp parallel for default(shared) private(j,id,dw)
    for (j = 0; j < size(); ++j) {
      dw = -1.0;
      id = omp_get_thread_num();

      if (vflag > 2) {
        cerr << setw(5) << j << " |> " << particle(j).getState() << " " 
             << setw(12) << "-" << " " 
             << setw(12) << particle(j).getWeight() << endl;
      }

      // cerr << "++ " << flush;
      dw = particle(j).force(eventTime,eventType,tree->ids[nextStep],(*rng)[id],pars);
//      if (dw == 0.0) {
//        const TransitionType* tt = model->getObsType(model->mapType(eventType));
//        char* str = new char[512];
//        sprintf(str,"Forced weight is zero is event '%s'.",tt->getName().c_str());
//        particle(j).msg(str);
//        delete[] str;
//      }
      // cerr << "++ " << flush;
      particle(j).updateWeight(dw);
      // cerr << "++ " << endl;

      if (vflag > 2) {
        cerr << setw(5) << j << "  > " << particle(j).getState() << " " 
             << setw(12) << dw << " " 
             << setw(12) << particle(j).getWeight() << endl;
      }

//      if (! noProb) {
//        /* probability that the event is included in the tree */
//        const TransitionType* tt = model->getObsType(modelType);
//        EpiState curState = particle(j).getState();
//        dw = tt->applyProb(curState,pars);
//        particle(j).updateWeight(dw);
//        particle(j).updateProb(dw);
//      }
    }

    return tree->ttypes[nextStep];
  }

  // ========================================================================

  void TrajParticleFilter::calcWeights(const void* pars) {
    size_t j;
#pragma omp parallel for default(shared) private(j)
    for (j = 0; j < size(); ++j) {
      double dw = particle(j).calcWeight(pars);
//      cerr << j << " " << particle(j).getProb() << " --> " << dw << endl;
      particle(j).updateProb(dw);
    }
  }

  // ========================================================================

  double TrajParticleFilter::meanWeight() const {
    double mean(0.0);
    for (size_t j(0); j < size(); ++j) mean += particle(j).getWeight();
    return mean/size();
  }

  // ========================================================================

  double TrajParticleFilter::meanProb() const {
    double mean(0.0);
    for (size_t j(0); j < size(); ++j) mean += particle(j).getProb();
    return log(mean/size());
  }

  // ========================================================================

  size_t TrajParticleFilter::stepTree(const Pars* pars, rng::Rng* rng, bool adjZero, double dt) 
  {
    // size_t nextStep = curStep + 1;
    size_t nextStep = curStep;
    if (nextStep < tree->times.size()) {
      double time = tree->times[nextStep];
      double step_time = -1.0;
      size_t j = 0;
      size_t id = 0;
      while (step_time < time) {
        step_time = cur_time() + dt;
        if (step_time > time) step_time = time;

#pragma omp parallel for default(shared) private(j,id)
        for (j = 0; j < size(); ++j) {
          id = omp_get_thread_num();
          particle(j).resetProb();
          int ret = particle(j).simulateTrajectory(step_time,pars,(*rng)[id]);
          if (ret < 0) {
            particle(j).setWeight(0.0);
          } else {
            particle(j).setWeight(particle(j).getProb());
          }
        }
        if (step_time < time) sampleInPlace((*rng)[0]);
      }
      return nextStep;
    } else {
      return -1;
    }
  }

  // ========================================================================

  void TrajParticleFilter::printNames(char sep) const {
    for (size_t j(0); j < size(); ++j) {
      cout << particle(j).getName() << " " << particle(j).getWeight() << sep;
    }
    cout << endl;
  }

  // ========================================================================

  void TrajParticleFilter::print(size_t j) const {;
    cout << "# " << particle(j).getName() << endl;
    cout << (Trajectory) particle(j) << endl;
  }

  // =========================================================================

  void TrajParticleFilter::printAll() const {
    for (size_t j(0); j < size(); ++j) print(j);
  }

  // =========================================================================

  int TrajParticleFilter::sampleInPlace(rng::RngStream* rng) {
    size_t n = size();
    size_t i;

    vector<double> w(n,0.0);
    double totalW = w_vec(w);

    if (vflag > 1) cout << "SIP :: Sum = " << totalW << endl;

    if (totalW <= 0.0) {
      logw = -INFINITY;
      return -1;
    }

    logw += log(totalW/n);

    vector<unsigned int> samples(n);
    vector<double> w0(n+1,0.0);
    rng::make_discrete(n,w.data(),w0.data());
    for (i = 0; i < n; ++i) samples[i] = rng->discrete_x(n,w0.data());

    vector<int> counts(n,0);
    for (i = 0; i < n; ++i) ++counts[samples[i]];

    set<int> free_pos;
    for (i = 0; i < n; ++i) {
      // check if the particle survives
      if (counts[i] > 0) {
        --counts[i];
      } else {
        free_pos.insert(free_pos.end(),i);
      }
    }

    i = 0;
    while (free_pos.size() > 0) {
      // find next replacement particle
      while (i < n && counts[i] <= 0) ++i;
      int j = *free_pos.begin();
      free_pos.erase(free_pos.begin());
      particle(j) = particle(i);
      --counts[i];
    }

    return i;
  }
}

