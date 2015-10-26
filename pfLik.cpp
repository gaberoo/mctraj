#include "pfLik.h"
#include "ascii.h"

namespace MCTraj {
  double pfLik(const Model* m, const EpiState& init, const Tree& tree,
               const PFPars& pars, rng::Rng* rng, Trajectory* out) 
  {
    if (pars.vflag > 1) cerr << "Starting..." << flush;

    Trajectory T(init,m);

    TrajParticleFilter* pf;
    if (pars.history) {
      if (pars.vflag > 1) cerr << "History filter." << endl;
      pf = new HistoryFilter(m);
    } else {
      if (pars.vflag > 1) cerr << "Static filter." << endl;
      pf = new StaticFilter(m);
    }
    pf->setVerbosity(pars.vflag);

    if (pars.vflag > 1) cerr << "done." << endl;

    // add particles
    for (size_t i = 0; i < pars.num_particles; ++i) {
      char name[32];
      sprintf(name,"P%lu",i);
      pf->push_back(TrajParticle(string(name),1.0,T));
      (*pf)[i].setId(i);
      (*pf)[i].storeTrans(false);
    }

    // pf->copyFromPrev();

    // assign tree
    //   note: skip doesn't have an effect for StaticFilter
    pf->setTree(&tree,pars.skip,rng);

    int ret = 1;
    int filterRet = 0;

    double time = 0.0;

    while (ret >= 0) {
      if (pars.vflag > 1) cerr << "Advancing tree..." << endl;
      ret = pf->stepTree(m->getPars(),rng,pars.adj_zero,pars.step_size);
      if (pars.vflag > 1) cerr << "Finished advancing tree." << endl;

      if (ret >= 0) {
//        if (vflag > 1) cerr << "Calculating weights..." << endl;
//        pf.calcWeights(m->getPars());

        if (pars.vflag > 1) cerr << "Adding tree event..." << endl;
        pf->addTreeEvent(m->getPars(),rng);
        if (pars.vflag > 1) cerr << "Finished adding tree event." << endl;

        // if (pars.print_particles) pf.printFromLast();

        pf->setLast();

        if (pars.vflag > 1) cerr << "Filter particles:" << endl;
        if (pf->cur_time() >= time + pars.filter_time) {
          filterRet = pf->filter(rng);
          time += pars.filter_time;
        } else {
          filterRet = pf->filter(rng,'c');
        }

        if (filterRet < 0) {
          if (pars.vflag > 0) {
            cerr << "\033[1;31m";
            cerr << pf->cur_time() << "/" << pf->maxTime() 
                 << ": Particle collapse => " << pf->getCurType()->getName() << " !";
            cerr << "\033[0m" << endl;
            // if (pars.vflag > 2) pf.printFromLast(m->getPars(),1);
          }
          return -INFINITY;
        }
      } 
//      else {
//        cerr << "End of tree!" << endl;
//      }
    }

    double log_lik = pf->est();

    if (tree.extant > 0) {
      if (pars.vflag) cerr << "There are " << tree.extant << " extant nodes." << endl;
      if (m->getRho() > 0.0) {
        // Sampling at present
        for (size_t i = 0; i < pars.num_particles; ++i) {
          double w = m->sample_rho((*pf)[i].getState(),(*rng)[i]);
          (*pf)[i].updateWeight(w);
        }

        filterRet = pf->filter(rng);

        if (filterRet < 0) {
          cerr << "\033[1;31m";
          cerr << "Something went wrong with the sampling step at the end.";
          cerr << "\033[0m" << endl;
          return -INFINITY;
        }

        log_lik = pf->est() /* + gsl_sf_lnfact(tree.extant) */;
      } else {
        cerr << "\033[1;31m";
        cerr << "Extant species but rho = 0 !";
        cerr << "\033[0m" << endl;
        return -INFINITY;
      }
    }

    // pf.printMeanTraj(cout);
    // pf.printFromFirst();
    // if (out != NULL) *out = pf.singleTraj((*rng)[0]);

//    {
//      vector<double> mu;
//      vector<double> s2;
//      pf.weights(mu,s2);
//      cerr << "MU  [" << mu.size() << "]: ";
//      for (size_t i = 0; i < mu.size(); ++i) cerr << setw(10) << mu[i] << " ";
//      cerr << endl;
//      cerr << "VAR [" << s2.size() << "]: ";
//      for (size_t i = 0; i < s2.size(); ++i) cerr << setw(10) << s2[i] << " ";
//      cerr << endl;
//    }

//    cerr << m->getPars()->to_json() << endl;
//    cerr << log_lik << endl;

    if (pars.vflag > 1) cerr << "Cleaning up..." << flush;
    delete pf;
    if (pars.vflag > 1) cerr << "done." << endl;

    return log_lik;
  }
}

