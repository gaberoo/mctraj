#include "TrajParticleFilter.h"

namespace MCTraj {
  void TrajParticleFilter::setTree(const Tree* tree, int skip, Rng* rng) { 
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

  // =========================================================================

  void TrajParticleFilter::inc() { 
    if (size() <= ++curStep) {
      pf.push_back(vector<TrajParticle>(size()));
    }
  }

  // =========================================================================

  int TrajParticleFilter::filter(gsl_rng** rng, char filter_type) 
  {
    size_t n(size());

    vector<double> w(n,0.0);
    double totalW = w_vec(w);
    if (vflag > 1) cerr << "FIL (" << filter_type << ") :: Sum = " << totalW << endl;

    if (totalW <= 0.0) {
      logw = -INFINITY;
      return -1;
    }

    // allocate space for new generation of particles
    if (vflag > 1) cerr << "Increasing array size..." << flush;
    inc();
    if (vflag > 1) cerr << "done." << endl;

    if (filter_type == 'c') {
      size_t i;
#pragma omp parallel for default(shared) private(i)
      for (i = 0; i < n; ++i) {
        particle(i).copy(pf[curStep-1][i]);
        particle(i).setId(i);
        particle(i).setParent(i);
      }
    } else {
      logw += log(totalW/n);

      size_t i;

      if (vflag > 1) cerr << "Preparing weight vector..." << flush;
      for (i = 0; i < n; ++i) {
        w[i] /= totalW;
        pf[curStep][i].setWeight(w[i]);
        if (i > 0) w[i] += w[i-1];
      }
      if (vflag > 1) cerr << "done." << endl;

#pragma omp parallel for default(shared) private(i)
      for (i = 0; i < n; ++i) {
        try {
          // sampling with replacement
          size_t id = omp_get_thread_num();
//          if (id == 0 && vflag > 1) cerr << "Getting random number..." << flush;
          double ran = gsl_rng_uniform(rng[id]);
//          if (id == 0 && vflag > 1) cerr << "done." << endl;
          int p = (int) (n*ran);
          if (ran > w[p]) {
            // find the first element that's larger than ran
            while (ran > w[p] && p < (int) n) ++p;
          } else {
            // find the first element that's smaller than ran
            while (ran <= w[p] && p >= 0) --p;
            ++p;
          }
          particle(i).copy(pf[curStep-1][p]);
          particle(i).setId(i);
          particle(i).setParent(p);
          particle(i).resetWeight();
          particle(i).setProb(1.0);
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
      // cerr << i << " " << w[i] << endl;
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

  TrajParticle TrajParticleFilter::sample(gsl_rng* rng) const {
    size_t n(size());
    vector<double> w(n);
    w_vec(w);
    gsl_ran_discrete_t* pp = gsl_ran_discrete_preproc(n,w.data());
    size_t j = gsl_ran_discrete(rng,pp);
    gsl_ran_discrete_free(pp);
    return particle(j);
  }

  // ========================================================================

  int TrajParticleFilter::addTreeEvent(const void* pars, Rng* rng, int noProb) {
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

    double* rand = (double*) malloc(size()*sizeof(double));
    rng->uniform(size(),rand,0.0,1.0);

#pragma omp parallel for default(shared) private(j,id,dw)
    for (j = 0; j < size(); ++j) {
      dw = -1.0;
      id = omp_get_thread_num();

      if (vflag > 2) {
        cerr << setw(5) << j << " |> " << particle(j).getState() << " " 
             << setw(12) << "-" << " " 
             << setw(12) << particle(j).getWeight() << endl;
      }

      dw = particle(j).force(eventTime,eventType,tree->ids[nextStep],rng[id],pars);
      particle(j).updateWeight(dw);

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
    free(rand);

    return tree->ttypes[nextStep];
  }

  // ========================================================================

  void TrajParticleFilter::setLast() {
    size_t j;
#pragma omp parallel for default(shared) private(j)
    for (j = 0; j < size(); ++j) {
      particle(j).setLast();
    }
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

  double TrajParticleFilter::meanCWeight() const {
    double mean(0.0);
    int cnt(0);
    for (size_t j(0); j < size(); ++j) {
      if (isfinite(particle(j).cw())) {
        mean += particle(j).cw();
        ++cnt;
      }
    }
    return mean/(double) cnt;
  }

  // ========================================================================

  void TrajParticleFilter::resetWeights() {
    for (size_t j(0); j < size(); ++j) {
      particle(j).resetWeight();
    }
  }

  // ========================================================================

  size_t TrajParticleFilter::stepTree(const void* pars, gsl_rng** rng, bool adjZero, double dt) 
  {
    // size_t nextStep = curStep + 1;
    size_t nextStep = curStep;
    if (nextStep < tree->times.size()) {
      double time = tree->times[nextStep];
      double step_time = -1.0;
      size_t j = 0;
      size_t id = 0;
      while (step_time < time) {
        step_time = particle(0).getTime() + dt;
        if (step_time > time) step_time = time;

#pragma omp parallel for default(shared) private(j,id)
        for (j = 0; j < size(); ++j) {
          id = omp_get_thread_num();
          particle(j).resetProb();
          int ret = particle(j).simulateTrajectory(step_time,pars,rng[id]);
          if (ret < 0) {
            particle(j).setWeight(0.0);
          } else {
            particle(j).setWeight(particle(j).getProb());
          }
        }
        if (step_time < time) sampleInPlace(rng[0]);
      }
      return nextStep;
    } else {
      return -1;
    }
  }

  // ========================================================================

  size_t TrajParticleFilter::stepAdd(const void* pars, gsl_rng** rng) {
    size_t j;
    size_t zeros = 0;
    size_t nextStep = curStep;
    if (nextStep < tree->times.size()) {
#pragma omp parallel for default(shared) private(j) reduction(+:zeros)
      for (j = 0; j < size(); ++j) {
        size_t id = omp_get_thread_num();
        int ret = 3;
        size_t lz = 0;
        while (ret > 0) {
          ret = stepAddTP(j,pars,rng[id]);
          cerr << j << " > " << ret << ": " << particle(j).getWeight() << endl;
          if (ret > 0) {
            ++lz;
            cerr << "Ilegal simulation: " << particle(j).getState() << endl;
            cerr << (Trajectory) particle(j) << endl;
            cerr << " :: Retrying (" << lz << ")..." << endl;
          } else break;
        }
        if (ret < 0) {
          cerr << "There was an error in step-adding particle '" << j << "'." << endl;
        }
        zeros += lz;
      }
      return zeros;
    } else {
      return -1;
    }
  }

  // ========================================================================

  int TrajParticleFilter::stepAddTP(size_t j, const void* pars, gsl_rng* rng) {
    // size_t nextStep = curStep + 1;
    size_t nextStep = curStep;
    if (nextStep < tree->times.size()) {
      double eventTime = tree->times[nextStep];

      // Simulate to next event
      int ret = particle(j).simulateTrajectory(eventTime,pars,rng);
      if (ret < 0) { particle(j).setWeight(0.0); return ret; }

      double dw = particle(j).calcWeight(pars);

      if (dw <= 0.0) return 1;  // simulated illegal event

      // Add observed event
      int eventType = tree->ttypes[nextStep];
      int modelType = model->mapType(eventType);
      if (modelType < 0) return -102;

      dw = particle(j).force(eventTime,eventType,tree->ids[nextStep],rng,pars);
      particle(j).updateWeight(dw);

      if (dw <= 0.0) return 2;  // forced illegal event

      return 0;
    } else {
      return -101;
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

  void TrajParticleFilter::printFromLast(const void* pars, int offset) const {
    ostringstream prefix;
    prefix << curStep;
    for (size_t j(0); j < size(); ++j) {
      particle(j).print_from_last(cout,prefix.str(),pars,offset);
    }
  }

  // =========================================================================

  void TrajParticleFilter::printFromFirst(const void* pars) const {
    size_t i;
    size_t j;
    for (i = 0; i < pf.size(); ++i) {
      ostringstream prefix;
      prefix << i;
      for (j = 0; j < pf[i].size(); ++j) {
        pf[i][j].print_from_first(cout,prefix.str(),pars);
        cout << endl;
      }
    }
  }

  // =========================================================================

  void TrajParticleFilter::meanState(double x[], size_t gen) const {
    size_t i, j;
    size_t n = model->n();
    EpiState y;
    for (i = 0; i < n; ++i) x[i] = 0.0;
    for (j = 0; j < size(); ++j) {
      y = pf[gen][j].getState();
      for (i = 0; i < n; ++i) x[i] += y[i];
    }
    for (i = 0; i < n; ++i) x[i] /= 1.0*size();
  }

  // =========================================================================

  void TrajParticleFilter::varState(double var[], const double mean[], size_t gen) const {
    size_t i, j;
    size_t n = model->n();
    EpiState y;
    double z;
    for (i = 0; i < n; ++i) var[i] = 0.0;
    for (j = 0; j < size(); ++j) {
      y = pf[gen][j].getState();
      for (i = 0; i < n; ++i) {
        z = y[i]-mean[i];
        var[i] += z*z;
      }
    }
    for (i = 0; i < n; ++i) var[i] /= 1.0*(size()-1);
  }

  // =========================================================================

  void TrajParticleFilter::printMeanTraj(ostream& out) const {
    // Output mean trajectory
    double* x = new double[size()];
    double* y = new double[size()];
    size_t i, j;
    for (i = 0; i < pf.size()-1; ++i) {
      meanState(x,i);
      varState(y,x,i);
      out << setw(12) << tree->times[i] << " ";
      for (j = 0; j < model->n(); ++j) {
        out << setw(8) << x[j] << " " << setw(8) << y[j] << " ";
      }
      out << endl;
    }
    delete[] x;
    delete[] y;
  }

  // =========================================================================

  Trajectory TrajParticleFilter::singleTraj(gsl_rng* rng) const {
    vector<const TrajParticle*> traj;
    traj.reserve(pf.size());
    vector< vector<TrajParticle> >::const_reverse_iterator t;
    int parent = gsl_rng_uniform_int(rng,pf.size());
    const TrajParticle* tp;
    for (t = pf.rbegin(); t != pf.rend(); ++t) {
      tp = &((*t)[parent]);
      traj.push_back(tp);
      parent = tp->getParent();
    }
    vector<const TrajParticle*>::const_reverse_iterator tpi = traj.rbegin();
    Trajectory single(**tpi);
    ++tpi;
    for (; tpi != traj.rend(); ++tpi) single += **tpi;
    return single;
  }

  // =========================================================================

  int TrajParticleFilter::sampleInPlace(gsl_rng* rng) {
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
    gsl_ran_multinomial(rng,n,n,w.data(),samples.data());

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
      particle(*(free_pos.begin())) = particle(i);
      --counts[i];
    }

    return i;
  }
}

