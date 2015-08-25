#include "TrajParticleFilter.h"

namespace MCTraj {
  void TrajParticleFilter::setTree(const Tree* tree, int skip, rng::Rng* rng) 
  { 
    this->tree = tree; 
  }

  // =========================================================================

  int TrajParticleFilter::filter(rng::Rng* rng, char filter_type) 
  {
    size_t n = size();

    // calculate total weight
    vector<double> w(n,0.0);
    double totalW = w_vec(w);

    if (vflag > 1) {
      cerr << ascii::yellow;

      // output verbose info on weight distribution
      double s2 = gsl_stats_variance(w.data(),1,w.size());
      size_t k = gsl_stats_max_index(w.data(),1,n);

      cerr << "  FIL (" << filter_type << ") :: Sum = " << totalW 
           << ", var(w) = " << sqrt(s2)/(totalW/n) 
           << endl;

      cerr << "  Best particle = " << particle(k).getState()
           << endl;
    }

    // check that the total weight is positive
    if (totalW <= 0.0) { logw = -INFINITY; return -1; }

    // allocate space for new generation of particles
    if (vflag > 1) cerr << "  Increasing array size..." << flush;
    inc();
    if (vflag > 1) cerr << "done." << endl;

    if (filter_type == 'c') {
// copy filter ===============================================================
      size_t i;
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i)
#endif
      for (i = 0; i < n; ++i) {
        particle(i).copy(prev(i));
        particle(i).setId(i);
        particle(i).setParent(i);
      }
    } else {
// bootstrap filter ==========================================================
      logw += log(totalW/n);

      size_t i;

      if (vflag > 1) cerr << "  Preparing weight vector..." << flush;
      vector<double> x(n,0.0);
      rng::make_discrete(n,w.data(),x.data());
      if (vflag > 1) cerr << "done." << endl;

      int p;

#if defined(_OPENMP)
      size_t id;
#pragma omp parallel for default(shared) private(i,id,p)
#endif
      for (i = 0; i < n; ++i) {
        try {
          // sampling with replacement
#if defined(_OPENMP)
          id = omp_get_thread_num();
          p = (*rng)[id]->discrete_x(n,x.data());
#else
          p = (*rng)[0]->discrete_x(n,x.data());
#endif
          // copyFromPrev(i,p);
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

    if (vflag > 1) cerr << ascii::end << "Done filtering." << endl;

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

  int TrajParticleFilter::addTreeEvent(const void* pars, rng::Rng* rng) 
  {
    // get information on the next event
    size_t nextStep = curStep;
    double eventTime = tree->times[nextStep];
    int eventType = tree->ttypes[nextStep];
    int modelType = model->mapType(eventType);

    // output some verbose information
    if (vflag > 1) {
      cerr << ascii::yellow;
      cerr << "Add tree event :: " << setw(12) << eventTime << " " 
           << "(" << eventType << "->" << modelType << ") :: ";
      for (size_t i(0); i < tree->ids[nextStep].size(); ++i) {
        cerr << tree->ids[nextStep][i] << " ";
      }
      cerr << ascii::end << endl;
    }

    // make sure the model is specified, otherwise stop
    if (modelType < 0) return tree->ttypes[nextStep];

    size_t j;
    double dw;
    int cnt_zero = 0;

#if defined(_OPENMP)
    int id;
#pragma omp parallel for default(shared) private(j,id,dw)
#endif
    for (j = 0; j < size(); ++j) {
      dw = -1.0;

      if (vflag > 2) {
        cerr << setw(5) << j << " |> " << particle(j).getState() << " " 
             << setw(12) << "-" << " " 
             << setw(12) << particle(j).getWeight() << endl;
      }

      // cerr << "++ " << flush;
      if (particle(j).getWeight() > 0.0) {
#if defined(_OPENMP)
        id = omp_get_thread_num();
        dw = particle(j).force(eventTime,eventType,tree->ids[nextStep],(*rng)[id],pars);
#else
        dw = particle(j).force(eventTime,eventType,tree->ids[nextStep],(*rng)[0],pars);
#endif
      }
#ifdef DEBUG
      else {
        cerr << ascii::red 
             << " [" << j << "] Weight already zero, not adding any event." 
             << ascii::end << endl;
      }
#endif

      if (dw == 0.0) {
#ifdef DEBUG
        cerr << ascii::red 
             << " [" << j << "] Weight is zero after adding event." 
             << ascii::end << endl;
#endif
        cnt_zero++;
      }

//      if (dw == 0.0) {
//        const TransitionType* tt = model->getObsType(model->mapType(eventType));
//        char* str = new char[512];
//        sprintf(str,"Forced weight is zero in event '%s'.",tt->getName().c_str());
//        particle(j).msg(str);
//        delete[] str;
//      }

      // update the particle weight, which at this point should be equal to
      // the (conditioned) trajectory weight
      particle(j).updateWeight(dw);

      // output verbose information
      if (vflag > 2) {
        cerr << setw(5) << j << "  > " << particle(j).getState() << " " 
             << setw(12) << dw << " " 
             << setw(12) << particle(j).getWeight() << endl;
      }
    }
    
    if (vflag > 0) {
      cerr << "Number of zero-weight particles = " << cnt_zero << endl;
    }

    return tree->ttypes[nextStep];
  }

  // ========================================================================

  void TrajParticleFilter::calcWeights(const void* pars) {
    size_t j;
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(j)
#endif
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

  size_t TrajParticleFilter::stepTree
    (const Pars* pars, rng::Rng* rng, bool adjZero, double dt) 
  {
    // get step number
    size_t nextStep = curStep;

    // make sure there are enough steps in the times file
    if (nextStep < tree->times.size()) 
    {
      double time = tree->times[nextStep];
      double step_time = -1.0;
      size_t j = 0;
      int ret = 0;
#if defined(_OPENMP)
      size_t id = 0;
#endif

      // iterate over all incremental steps
      while (step_time < time) {
        step_time = cur_time() + dt;
        if (step_time > time) step_time = time;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(j,id,ret)
#endif
        for (j = 0; j < size(); ++j) {
          // set probability of the trajectory to zero
          // particle(j).resetProb();

          // simulate a trajectory
#if defined(_OPENMP)
          id = omp_get_thread_num();
          ret = particle(j).simulateTrajectory(step_time,pars,(*rng)[id],adjZero);
#else
          ret = particle(j).simulateTrajectory(step_time,pars,(*rng)[0],adjZero);
#endif

#ifdef DEBUG
          cerr << ascii::green
               << "  Particle (" << particle(j).getId() << "/" 
               << particle(j).getParent() << "): " 
               << particle(j).getState() 
               << ascii::end << endl;
#endif

          // check if the simulation was successful
          if (ret < 0) {
            // if no, set particle weight to zero
            particle(j).setWeight(0.0);
          } else {
            // if yes, get the probabiliy of the trajectory,
            // potentially conditioned on some illegal events
            // not happening
            particle(j).setWeight(particle(j).getProb());
          }
        }

        // if there are incremental steps, resample the particles in place
        // if (step_time < time) sampleInPlace((*rng)[0]);
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
    vector<double> w0(n,0.0);
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

