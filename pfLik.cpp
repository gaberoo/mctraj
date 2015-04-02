#include "pfLik.h"

namespace MCTraj {
  double pfLik(const Model* m, const EpiState& init, const Tree& tree,
               const PFPars& pars, rng::Rng* rng, Trajectory* out) 
  {
    Trajectory T(init,m);
    TrajParticleFilter pf(m);
    pf.setVerbosity(pars.vflag);

    for (size_t i(0); i < pars.num_particles; ++i) {
      char name[32];
      sprintf(name,"P%lu",i);
      // cerr << name << endl;
      pf.push_back(TrajParticle(string(name),1.0,T));
      pf[i].setId(i);
    }
    pf.setTree(&tree,pars.skip,rng);

    int ret = 1;
    int filterRet = 0;

    double time = 0.0;

    while (ret >= 0) {
      if (pars.vflag > 1) cerr << "Advancing tree..." << endl;
//      ret = pf.stepAdd(m->getPars(),rng);
      ret = pf.stepTree(m->getPars(),rng);
      if (ret >= 0) {
//        if (vflag > 1) cerr << "Calculating weights..." << endl;
//        pf.calcWeights(m->getPars());

        if (pars.vflag > 1) cerr << "Adding tree event..." << endl;
        pf.addTreeEvent(m->getPars(),rng,1);

        if (pars.print_particles) pf.printFromLast();

        pf.setLast();

        if (pf[0].getTime() >= time + pars.filter_time) {
          if (pars.vflag > 1) cerr << "Filter particles..." << endl;
          filterRet = pf.filter(rng);
          time += pars.filter_time;
        } else {
          filterRet = pf.filter(rng,'c');
        }

        if (filterRet < 0) {
          if (pars.vflag > 0) {
            cerr << "\033[1;31m";
            cerr << pf[0].getTime() << "/" << pf.maxTime() 
                 << ": Particle collapse !";
            cerr << "\033[0m" << endl;
            if (pars.vflag > 2) pf.printFromLast(m->getPars(),1);
          }
          return -INFINITY;
        }
      }
    }

    double log_lik = pf.est();

    if (tree.extant > 0) {
      if (m->getRho() > 0.0) {
        // Sampling at present
        for (size_t i = 0; i < pars.num_particles; ++i) {
          double w = m->sample_rho(pf[i].getState(),rng[i]);
          pf[i].updateWeight(w);
        }

        filterRet = pf.filter(rng);

        if (filterRet < 0) {
          cerr << "\033[1;31m";
          cerr << "Something went wrong with the sampling step at the end.";
          cerr << "\033[0m" << endl;
          return -INFINITY;
        }

        log_lik = pf.est() /* + gsl_sf_lnfact(tree.extant) */;
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

