#include "HistoryFilter.h"

namespace MCTraj {
  void MCTraj::HistoryFilter::setTree(const Tree* tree, int skip, rng::Rng* rng) {
    TrajParticleFilter::setTree(tree,skip,rng);

    // preallocate the space for the particles
    if (vflag > 1) cerr << "allocating..." << flush;
    pf.reserve(tree->size()+1);
    for (size_t i = 0; i <= tree->size(); ++i) {
      pf.push_back(vector<TrajParticle>(size()));
    }
    if (vflag > 1) cerr << "done." << endl;

    // add tree events
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

  void HistoryFilter::printMeanTraj(ostream& out) const {
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

  Trajectory HistoryFilter::singleTraj(rng::RngStream* rng) const {
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

  void HistoryFilter::printFromLast(const void* pars, int offset) const {
    ostringstream prefix;
    prefix << curStep;
    for (size_t j(0); j < size(); ++j) {
      particle(j).print_from_last(cout,prefix.str(),pars,offset);
    }
  }

  // =========================================================================

  void HistoryFilter::printFromFirst(const void* pars) const {
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

  void HistoryFilter::meanState(double x[], size_t gen) const {
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

  void HistoryFilter::varState(double var[], const double mean[], size_t gen) const {
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

  // ========================================================================

  void HistoryFilter::weights(vector<double>& mu, vector<double>& s2) const {
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

  // =========================================================================

  void HistoryFilter::inc() { 
    if (size() <= ++curStep) {
      pf.push_back(vector<TrajParticle>(size()));
    }
  }

  // ========================================================================

  void HistoryFilter::setLast() {
    size_t j;
#pragma omp parallel for default(shared) private(j)
    for (j = 0; j < size(); ++j) {
      particle(j).setLast();
    }
  }

}

