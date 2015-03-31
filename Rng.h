#ifndef __RNG_H__
#define __RNG_H__

#ifdef MKLRNG
#include "mkl_vsl.h"
#endif

#include <gsl/gsl_rng.h>
#include <omp.h>

namespace MCTraj {
  class Rng {
    public:
      Rng() {}
      virtual ~Rng() {}

      virtual void alloc() = 0;
      virtual void free() = 0;
      virtual void uniform(size_t n, double* r, double a, double b) = 0;
  }

#ifdef MKLRNG
  class MKLRng : public Rng {
    public:
      MKLRng() {}
      virtual MKLRng() {}

      inline void alloc() {
        vslNewStream(&stream,VSL_BRNG_MT19937,777);
      }

      inline void free() {
        vslDeleteStream(&stream); 
      }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,r,a,b);
      }

    protected:
      VSLStreamStatePtr stream;
  };
#endif

  class GSLRng : public Rng {
    public:
      GSLRng() : size(0), rng(NULL) {}
      virtual ~GSLRng() { free(); }

      inline void alloc(size_t n, const gsl_rng_type* T = gsl_rng_taus2) {
        rng = (gsl_rng**) malloc(n*sizeof(gsl_rng*));
        for (size_t i = 0; i < n; ++i) {
          rng[i] = gsl_rng_alloc(gsl_rng_taus2);
          gsl_rng_set(rng[i],seed+i);
        }
        size = n;
      }

      inline void free() {
        if (size > 0) {
          for (size_t i = 0; i < size; ++i) gsl_rng_free(rng[i]);
          free(rng);
        }
      }

      inline gsl_rng* operator[](size_t i) { return rng[i]; }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        size_t i;
#pragma omp parallel for default(shared) private(i)
        for (i = 0; i < n; ++i) {
          size_t id = omp_get_thread_num();
          if (id < size) {
            r[i] = a + (b-a)*gsl_rng_uniform(rng[id]);
          } else {
            r[i] = 0.0;
          }
        }
      }

    protected:
      size_t size;
      gsl_rng** rng;
  };
}

#endif

