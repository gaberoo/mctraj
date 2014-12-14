#include "pfLik.h"

namespace MCTraj {
  double pfLik(const Model* m, const EpiState& init, const Tree& tree,
               size_t num_particles, gsl_rng** rng, int vflag,
               Trajectory* out, int skip, int print_particles) 
  {
    Trajectory T(init,m);
    TrajParticleFilter pf(m);
    pf.setVerbosity(vflag);

    for (size_t i(0); i < num_particles; ++i) {
      char name[32];
      sprintf(name,"P%lu",i);
      // cerr << name << endl;
      pf.push_back(TrajParticle(string(name),1.0,T));
      pf[i].setId(i);
    }
    pf.setTree(&tree,skip,rng);

    int ret = 1;
    int filterRet = 0;

    while (ret >= 0) {
      if (vflag > 1) cerr << "Advancing tree..." << endl;
      ret = pf.stepTree(m->getPars(),rng);
      if (ret >= 0) {
        if (vflag > 1) cerr << "Calculating weights..." << endl;
        pf.calcWeights(m->getPars());

        if (vflag > 1) cerr << "Adding tree event..." << endl;
        pf.addTreeEvent(m->getPars(),rng,1);

        if (print_particles) pf.printFromLast();

        pf.setLast();

        if (vflag > 1) cerr << "Filter particles..." << endl;
        filterRet = pf.filter(rng);

        if (filterRet < 0) {
          if (vflag > 0) {
            cerr << "\033[1;31m";
            cerr << pf[0].getTime() << "/" << pf.maxTime() 
                 << ": Particle collapse !";
            cerr << "\033[0m" << endl;
            if (vflag > 1) pf.printFromLast(m->getPars(),1);
          }
          return -INFINITY;
        }
      }
    }

    double log_lik = pf.est();

    if (tree.extant > 0) {
      if (m->getRho() > 0.0) {
        // Sampling at present
        for (size_t i = 0; i < num_particles; ++i) {
          size_t k = pf[i].getState(2);
          size_t I = pf[i].getState(1);
          double w = 1.0;
          if (m->getRho() >= 1.0 && I == k) {
            w = 1.0;
          } else if (m->getRho() < 1.0 && I >= k) {
            // w = gsl_pow_int(pars.rho,k) * gsl_pow_int(1-pars.rho,I-k);
            w = gsl_ran_binomial_pdf(k,m->getRho(),I);
          } else {
            w = 0.0;
          }
          // cout << i << " " << w << endl;
          pf[i].updateWeight(w);
        }

        filterRet = pf.filter(rng);
        if (filterRet < 0) {
          cerr << "\033[1;31m";
          cerr << "Something went wrong with the sampling step at the end.";
          cerr << "\033[0m" << endl;
          return -INFINITY;
        }

        log_lik = pf.est() + gsl_sf_lnfact(tree.extant);
     } else {
        cerr << "\033[1;31m";
        cerr << "Extant species but rho = 0 !";
        cerr << "\033[0m" << endl;
        return -INFINITY;
      }
    }

    // pf.printMeanTraj(cout);
    // pf.printFromFirst();
    if (out != NULL) *out = pf.singleTraj(rng[0]);

    return log_lik;
  }
}

