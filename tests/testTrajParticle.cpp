#include "TrajParticleFilter.h"
using namespace MCTraj;

#include <iostream>
#include <string>
#include <cstring>
using namespace std;

struct EpiPars {
  double N;
  double beta;
  double gamma;
  double psi;
};

/****************************************************************************/

double infRateFun(const EpiState& es, void* pars) {
  EpiPars ep = *(EpiPars*) pars;
  return ep.beta*es[0]*es[1]/ep.N;
}

double recovRateFun(const EpiState& es, void* pars) {
  EpiPars ep = *(EpiPars*) pars;
  return (ep.gamma+ep.psi)*es[1];
}

/****************************************************************************/

double treeProbInf(const EpiState& es, const vector<int>& nlin, void* pars) 
{
  double w = 1.0;
  int ninf = es[1];
  if (ninf >= nlin[0]) {
    w *= 1.*nlin[0]*(nlin[0]-1.0)/(ninf*(ninf-1.0));
  } else {
    w = 0.0;
  }
  return w;
}

double treeProbRecov(const EpiState& es, const vector<int>& nlin, void* pars) 
{
  EpiPars ep = *(EpiPars*) pars;
  double s = ep.psi/(ep.psi+ep.gamma);
  double I = es[1]+1.0;
  double k = nlin[0]+1.0;
  return (I >= k) ? s*k/((1.-s)*(I-k)+s*k) : 0.0;
}

/****************************************************************************/

void treeTransFromType(int eventType, vector<int>& x) {
  if (x.size() < 1) x.resize(1);
  switch (eventType) {
    case 1:
      x[0] = 1;
      break;

    case 0:
      x[0] = -1;
      break;

    default:
      x[0] = 0;
      break;
  }
}

/****************************************************************************/

int main(int argc, char** argv) {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,time(NULL));

  const size_t nstates = 2;
  const int infChange[] = { -1, 1 };
  const int recoverChange[] = { 1, -1 };

  EpiPars pars = { 100.0, 1.0, 0.1, 0.1 };
  size_t num_particles = 1;
  int reps = 1;

  int vflag = 0;
  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"N:b:u:s:vn:r:")) != -1) {
    switch (c) {
      case 'N': pars.N = atof(optarg); break;
      case 'b': pars.beta = atof(optarg); break;
      case 'u': pars.gamma = atof(optarg); break;
      case 's': pars.psi = atof(optarg); break;
      case 'v': ++vflag; break;
      case 'n': num_particles = atoi(optarg); break;
      case 'r': reps = atoi(optarg); break;
    }
  }

  EpiState es(nstates);
  es[1] = 1;
  es[0] = (int) (pars.N - es[1]);

  Tree tree(*(argv+optind));
  tree.reverse();

  Trajectory T(es);

  T.addTransType(nstates,recoverChange,recovRateFun,treeProbRecov);
  T.addTransType(nstates,infChange,infRateFun,treeProbInf);

  for (int r(0); r < reps; ++r) {
    TrajParticleFilter pf;
    pf.set_etype_fun(treeTransFromType);

    for (size_t i(0); i < num_particles; ++i) {
      char name[32];
      sprintf(name,"P%d",(int) i);
      pf.push_back(TrajParticle(string(name),1.0,T));
    }
    pf.setTree(&tree);
    pf.inc();

    int ret = 1;
    int filterRet;

    while (ret > 0) {
      ret = pf.stepTree(&pars,rng);
      if (ret > 0) {
        pf.calcWeights(&pars);
        pf.printFromLast();
        pf.addTreeEvent(&pars);
        // pf.observeTreeEvent(INFINITY,&pars);
        filterRet = pf.filter(rng);
        pf.resetWeights();
        if (filterRet < 0) {
          cerr << pf[0].getTime() << "/" << pf.maxTime() 
               << ": " << "Particle collapse!" << endl;
          break;
        }
//        for (size_t i = 0; i < pf.size(); ++i) {
//          cout << i << " ";
//          pf[i].print_from_last(cout);
//        }
//        cout << endl << endl;
      }
      // pf.inc();
    }

    cerr << "LL(w) = " << pf.est() << endl;
    cerr << "LL(p) = " << pf.est_p() << endl;

    vector<int> init(1,1);
    cerr << pf[0].cw() << " " 
         << log(pf.trajProb(0,0.0,INFINITY,init,&pars)) << endl;

    // pf.print(gsl_rng_uniform_int(rng,pf.size()));
  }

  gsl_rng_free(rng);
  return 0;
}

