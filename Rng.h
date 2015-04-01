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
      Rng() : isAlloc(false), seed(time(NULL)) {}
      virtual ~Rng() {}

      virtual void alloc() = 0;
      virtual void free() = 0;
      virtual void uniform(size_t n, double* r, double a, double b) = 0;

      inline void set_seed(unsigned long s) { seed = s; }

    protected:
      bool isAlloc;
      unsigned long seed;
  };

#ifdef MKLRNG
  class MKLRng : public Rng {
    public:
      MKLRng() {}
      virtual ~MKLRng() { free(); }

      inline void alloc() {
        if (! isAlloc) {
          vslNewStream(&stream,VSL_BRNG_MT19937,seed);
          isAlloc = true;
        }
      }

      inline void free() {
        if (isAlloc) {
          vslDeleteStream(&stream); 
          isAlloc = false;
        }
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
      GSLRng(size_t n = 1) : size(n), type(gsl_rng_taus2), rng(NULL) {}

      virtual ~GSLRng() { this->free(); }

      inline void alloc() {
        if (! isAlloc) {
          rng = (gsl_rng**) malloc(size*sizeof(gsl_rng*));
          for (size_t i = 0; i < size; ++i) {
            rng[i] = gsl_rng_alloc(gsl_rng_taus2);
          }
          isAlloc = true;
          rng_seed();
        }
      }

      inline void free() {
        if (isAlloc) {
          for (size_t i = 0; i < size; ++i) gsl_rng_free(rng[i]);
          std::free(rng);
          rng = NULL;
          isAlloc = false;
        }
      }

      inline gsl_rng* operator[](size_t i) { return rng[i]; }
      inline gsl_rng** ptr() { return rng; }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        size_t i;
#pragma omp parallel for default(shared) private(i) num_threads(size)
        for (i = 0; i < n; ++i) {
          size_t id = omp_get_thread_num();
          if (id < size) {
            r[i] = a + (b-a)*gsl_rng_uniform(rng[id]);
          } else {
            cerr << "RNG array not large enough: " << id << " of " << size << "." << endl;
            r[i] = 0.0;
          }
        }
      }

      inline void set_size(size_t n) { size = n; }
      inline void set_type(const gsl_rng_type* T) { type = T; }
      inline void rng_seed() {
        for (size_t i = 0; i < size; ++i) {
          gsl_rng_set(rng[i],seed+i);
        }
      }

    protected:
      size_t size;
      const gsl_rng_type* type;
      gsl_rng** rng;
  };
}

#endif

