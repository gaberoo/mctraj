#ifndef __RNGSTREAM_H__
#define __RNGSTREAM_H__

namespace Rng {
  class RngStream {
    public:
      RngStream() : is_alloc(0) {}
      virtual ~RngStream() {}
      virtual void alloc(unsigned long seed) = 0;
      virtual void free() = 0;
      virtual void uniform(size_t n, double* r, double a, double b) = 0;

      inline void leapfrog(size_t k, size_t n) {}

    protected:
      int is_alloc;
  };
}

#endif // __RNGSTREAM_H__

