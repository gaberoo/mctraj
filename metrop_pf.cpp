#include <stdlib.h>
#include <iostream>
#include <vector>
#include <ctime>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "../Forest.h"
#include "pfLik.h"

int maxit = 1e6;
int nvar = 3;
int nstates = 2;
size_t num_particles = 100;
gsl_rng* rng = NULL;

vector<double> x(3);

double priorFun(const vector<double>& y, const vector<double>& s) {
  double p = 0.0;
  p += log(gsl_ran_lognormal_pdf(y[0],log(x[0]),s[0]));
  p += log(gsl_ran_gaussian_pdf(y[1]-x[1],s[1]));
  p += log(gsl_ran_gaussian_pdf(y[2]-x[2],s[2]));
  return p;
}

double likFun(const vector<double>& y, const EpiState& init, const Tree& tree) {
  EpiPars pars = { y[0], y[1], y[2], y[2] };
  return pfLik_SIS(pars,init,tree,num_particles,rng);
}

int main(int argc, char* argv[]) {
  // if (parseArgs(argc,argv) == -1) return 0;
  if (optind == argc) {
    cerr << "Please provide at least one tree file.\n" << endl;
    return 0;
  }

  optind = 1;
  Tree tree(*(argv+optind));
  tree.reverse();

  // =========================================================================

  rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng,time(NULL));

  EpiState es(2);

  double last_lik = -INFINITY;
  double new_lik = last_lik;

  double last_prior = -INFINITY;
  double new_prior = last_prior;

  x[0] = 100.0;
  x[1] = 1.0;
  x[2] = 0.1;

  vector<double> sigma(x);
  vector<double> last(x);
  last[0] = gsl_ran_gaussian(rng,sqrt(sigma[0]))+last[0];
  last[1] = gsl_ran_gaussian(rng,sqrt(sigma[1]))+last[1];
  last[2] = gsl_ran_gaussian(rng,sqrt(sigma[2]))+last[2];

  last_prior = priorFun(last,sigma);
  cerr << "Prior     = " << last_prior << endl;

  es[1] = 1;
  es[0] = (int) (last[0] - es[1]);
  last_lik = likFun(last,es,tree);

  cerr << "Liklihood = " << last_lik   << endl;

  if (! gsl_isinf(last_lik)) {
    for (int it = 1; it < maxit; ++it) {
      for (int j = 0; j < nvar; ++j) {
        int accept = 0;
        vector<double> proposal(last);
        proposal[j] = gsl_ran_gaussian(rng,sqrt(sigma[j])) + last[j];
        new_prior = priorFun(proposal,sigma);
        if (! gsl_isinf(new_prior)) {
          es[1] = 1;
          es[0] = (int) (proposal[0] - es[1]);
          new_lik = likFun(proposal,es,tree);
          double alpha = exp(new_lik + new_prior - last_lik - last_prior);
          if (alpha > gsl_rng_uniform(rng)) accept = 1;
          if (accept) {
            last = proposal;
            last_lik = new_lik;
            last_prior = new_prior;
          }
        }
      }
      printf("%f %f %f %f %f\n",last[0],last[1],last[2],last_lik,last_prior);
    }
  }

  gsl_rng_free(rng);
  return 0;
}


