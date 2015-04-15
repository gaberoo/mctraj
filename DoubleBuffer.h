#ifndef __DOUBLEBUFFER_H__
#define __DOUBLEBUFFER_H__

template<typename T>
class DoubleBuffer {
  protected:
    size_t n;
    T* first;
    T* second;
    T* data;

  public:
    DoubleBuffer(size_t n) { alloc(n); }
    virtual ~DoubleBuffer() { free(); }

    inline void swap() { T* tmp = first; first = second; second = tmp; }
    inline T* one() { return first; }
    inline T* two() { return second; }

  protected:
    inline void alloc(size_t n) { 
      free();
      data = new T[2*n];
      first = data;
      second = data + n;
      this->n = n;
    }

    inline void free() {
      if (n > 0) {
        delete[] data;
        data = NULL;
        first = NULL;
        second = NULL;
        n = 0;
      }
    }
};

#endif
