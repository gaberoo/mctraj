#ifndef __MKLSTREAM_H__
#define __MKLSTREAM_H__

#include <mkl_vsl.h>
#include "RngStream.h"

namespace Rng {
  class MKLStream : public RngStream {
    public:
      MKLStream() : brng(VSL_BRNG_MCG31) {}
      virtual ~MKLStream() { free(); }

      inline void alloc(unsigned long seed = time(NULL)) {
        if (! is_alloc) {
          vslNewStream(&stream,brng,seed);
          is_alloc = 1;
        }
      }

      inline void free() {
        if (is_alloc) {
          vslDeleteStream(&stream); 
          is_alloc = 0;
        }
      }

      inline void copy(const MKLStream& s) {
        if (is_alloc) free();
        vslCopyStream(&stream,s.stream);
        is_alloc = 1;
      }

      inline void leapfrog(const MKL_INT k, const MKL_INT n) {
        vslLeapfrogStream(stream,k,n); 
      }

      inline void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,r,a,b);
      }

      inline void set_type(const MKL_INT type) { brng = type; }

    protected:
      MKL_INT brng;
      VSLStreamStatePtr stream;
  };
}

#endif

