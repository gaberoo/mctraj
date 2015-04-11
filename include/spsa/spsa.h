#ifndef __SPSA_H__
#define __SPSA_H__

#include <gsl/gsl_math.h>
#include <rng/RngStream.h>

namespace spsa {
  typedef double (*LossFun)(const double* state, void* pars);

  typedef struct {
    size_t n;
    double a;
    double c;
    double alpha;
    double gamma;
    double A;
    double ak;
    double ck;
    double* y;
    double* x1;
    double* x2;
    double* delta;
    double* grad;
    double* lo;
    double* hi;
    rng::RngStream* rng;
    LossFun fun;
    void* pars;
  } pars_t;

  inline void init_pars(pars_t* p, size_t n) {
    p->n = n;
    p->y     = new double[2];
    p->lo    = new double[n];
    p->hi    = new double[n];
    p->grad  = new double[n];
    p->x1    = new double[n];
    p->x2    = new double[n];
    p->delta = new double[n];
  }

  inline void free_pars(pars_t* p) {
    delete[] p->delta;
    delete[] p->x2;
    delete[] p->x1;
    delete[] p->y;
    delete[] p->grad;
    delete[] p->lo;
    delete[] p->hi;
  }

  inline int approx_gradient(const pars_t* p, const double* theta) {
    // cerr << "Approximating gradient..." << endl;
    double r = 0.0;

    for (size_t i = 0; i < p->n; ++i) {
      p->rng->uniform(1,&r);
      p->delta[i] = (r >= 0.5) ? 1.0 : -1.0;
      p->x1[i] = theta[i] + p->ck*p->delta[i];
      p->x2[i] = theta[i] - p->ck*p->delta[i];
    }

    // cerr << "  ... calculating functions ..." << endl;
    p->y[0] = p->fun(p->x1,p->pars);
    p->y[1] = p->fun(p->x2,p->pars);

    int f1 = gsl_finite(p->y[1]);
    int f2 = gsl_finite(p->y[2]);

    if (! (f1 & f2)) return -1;

    int ret = 3;

    if (f1 ^ f2) {
      // one-sided approximation
      // cerr << " ... one-sided approximation ..." << endl;
      if (f1) {
        p->y[1] = p->fun(theta,p->pars);
        ret = 1;
      } else {
        p->y[0] = p->fun(theta,p->pars);
        ret = 2;
      }
      for (size_t i = 0; i < p->n; ++i) p->delta[i] *= 0.5;
    }

    double dy = p->y[0]-p->y[1];
    for (size_t i = 0; i < p->n; ++i) {
      p->grad[i] = dy/(2.*p->ck*p->delta[i]);
    }

    return ret;
  }

//  inline int run_spsa(int steps, size_t n, double* theta) {
//    double ak = 0.0;
//    double ck = 0.0;
//    double* grad = new double[n];
//
//    size_t i;
//    size_t k = 0;
//    // printf("%d %f\n",k,th[0]);
//
//    for (k = 1; k <= steps; ++k) {
//      ak = p.a/pow(k+1.0+p.A,p.alpha);
//      ck = p.c/pow(k+1.0,p.gamma);
//      approx_gradient(1,th,ck,grad,&p);
//      for (i = 0; i < n; ++) theta[i] -= ak*grad[i];
//    }
//
//    delete[] grad;
//  }
}

#endif 
