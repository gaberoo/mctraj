#include "TrajParticleFilter.h"

namespace MCTraj {
  void TrajParticleFilter::setTree(const Tree* tree, int skip, rng::Rng* rng) { 
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

  int TrajParticleFilter::filter(rng::Rng* rng, char filter_type) 
  {
    size_t n(size());

    vector<double> w(n,0.0);
    double totalW = w_vec(w);
    if (vflag > 1) {
      double s2 = gsl_stats_variance(w.data(),1,w.size());
      cerr << "FIL (" << filter_type << ") :: Sum = " << totalW 
           << ", var(w) = " << sqrt(s2)/(totalW/n) << endl;
      size_t k = gsl_stats_max_index(w.data(),1,n);
      cerr << "Best particle = " << particle(k).getState() << endl;
    }

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
          particle(i).copy(pf[curStep-1][p]);
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

  void TrajParticleFilter::weights(vector<double>& mu, vector<double>& s2) const {
    mu.clear();
    s2.clear();
    mu.reserve(pf.size());
    s2.reserve(pf.size());
    double wtot;
    size_t n = pf[0].size();
    vector<double> w(n);
    for (size_t m = 0; m < pf.size() && m <= curStep; ++m) {
      n = pf[m].size();
      w.resize(n);
      wtot = 0.0;
      for (size_t i = 0; i < n; ++i) {
        w[i] = pf[m][i].getWeight();
        wtot += w[i];
      }
      mu.push_back(wtot/n);
      s2.push_back(gsl_stats_variance_m(w.data(),1,n,mu.back()));
    }
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
        step_time = particle(0).getTime() + dt;
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

  size_t TrajParticleFilter::stepAdd(const void* pars, rng::Rng* rng) {
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
          ret = stepAddTP(j,pars,(*rng)[id]);
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

  int TrajParticleFilter::stepAddTP(size_t j, const void* pars, rng::RngStream* rng) {
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

  Trajectory TrajParticleFilter::singleTraj(rng::RngStream* rng) const {
    vector<size_t> seq(pf.size(),0);
    int p = 0;
    size_t j = curStep;
    rng->uniform_int(1,&p,0,pf[j].size());
    seq[j] = p;
    while (j > 0) {
      seq[j-1] = pf[j][p].getParent();
      --j;
    }
    Trajectory out(pf[0][seq[0]]);
    for (j = 1; j <= curStep; ++j) {
      out += pf[j][seq[j]];
    }
    return out;
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

