#ifndef __SPSA_H__
#define __SPSA_H__

#include <rng/RngStream.h>

namespace spsa {
  typedef double (*LossFun)(double* state, void* pars);

  typedef struct {
    size_t n;
    double a;
    double c;
    double alpha;
    double gamma;
    double A;
    double ak;
    double ck;
    double* x1;
    double* x2;
    double* delta;
    double* grad;
    rng::RngStream* rng;
    LossFun fun;
    void* pars;
  } pars_t;

  inline int approx_gradient(const pars_t* p, const double* theta, void* pars) {
    double r = 0.0;

    for (size_t i = 0; i < n; ++i) {
      p->rng->uniform(1,&r);
      p->delta[i] = (r >= 0.5) ? 1.0 : -1.0;
      p->x1[i] = theta[i] + p->ck*p->delta[i];
      p->x2[i] = theta[i] - p->ck*p->delta[i];
    }

    double y1 = p->fun(p->x1,p->pars);
    double y2 = p->fun(p->x2,p->pars);
    double dy = y1-y2;

    for (size_t i = 0; i < n; ++i) {
      p->grad[i] = dy/(2.*p->ck*p->delta[i]);
    }

    return 0;
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
