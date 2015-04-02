#ifndef __RNGSTREAM_H__
#define __RNGSTREAM_H__

#include <cstdlib>
#include <ctime>

namespace rng {
  inline void make_discrete(size_t n, const double* w, double* x) {
    x[0] = 0.0;
    for (size_t i = 1; i < n; ++i) x[i] = x[i-1] + w[i-1];
    x[n] = x[n-1]+w[n-1];
  }

  class RngStream {
    public:
      RngStream() : is_alloc(0) {}
      virtual ~RngStream() {}
      virtual void alloc(unsigned long seed = time(NULL)) = 0;
      virtual void free() = 0;

      inline void leapfrog(size_t k, size_t n) {}

      virtual void uniform(size_t n, double* r, double a = 0.0, double b = 1.0) = 0;
      virtual void uniform_int(size_t n, int* r, int a = 0, int b = 100) = 0;

      int discrete(size_t n, const double* w) {
        double* w0 = new double[n+1];
        make_discrete(n,w,w0);
        size_t j = discrete_x(n,w0);
        delete[] w0;
        return j;
      }

      int discrete_x(size_t n, const double* w) {
        double r;
        uniform(1,&r,0,w[n]);
        size_t j = r/w[n]*n;
        while (j >= 0 && j < n) {
          if (r < w[j]) --j;
          else if (r >= w[j+1]) ++j;
          else break;
        }
        return j;
      }

    protected:
      int is_alloc;
  };
}

#endif // __RNGSTREAM_H__

