#ifndef __GSLSTREAM_H__
#define __GSLSTREAM_H__

#include "RngStream.h"
#include <gsl/gsl_rng.h>

namespace rng {
  class GSLStream : public RngStream {
    public:
      GSLStream() : type(gsl_rng_taus2) {}
      virtual ~GSLStream() { free(); }

      inline void alloc(unsigned long seed = time(NULL)) {
        if (! is_alloc) {
          rng = gsl_rng_alloc(type);
          is_alloc = 1;
          gsl_rng_set(rng,seed);
        }
      }

      inline void free() {
        if (is_alloc) {
          gsl_rng_free(rng);
          is_alloc = 0;
        }
      }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        for (size_t i = 0; i < n; ++i) {
          r[i] = a + (b-a)*gsl_rng_uniform(rng);
        }
      }

      inline void uniform_int(size_t n, int* r, int a = 0, int b = 10) {
        for (size_t i = 0; i < n; ++i) {
          r[i] = a + gsl_rng_uniform_int(rng,b-a);
        }
      }


    protected:
      const gsl_rng_type* type;
      gsl_rng* rng;
  };
}

#endif

