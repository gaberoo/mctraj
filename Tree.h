#ifndef __TREE_H__
#define __TREE_H__

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

#include <gsl/gsl_sort.h>

class Tree {
public:
  /* CONSTRUCTORS-DESTRUCTORS */

  // 1. Default
  Tree() : extant(0), maxExtant(0), nroot(0), fn(""), is_rev(0), numBranches(0) {}

  // 2. Read tree from file
  Tree(const char* fn) : nroot(0), is_rev(0) { read(string(fn)); }
  Tree(string fn) : nroot(0), is_rev(0) { read(fn); }

  // 3. Copy constructor
  Tree(const Tree& T) 
    : times(T.times), ttypes(T.ttypes), extant(T.extant),
      maxExtant(T.maxExtant), nroot(T.nroot), is_rev(T.is_rev),
      numBranches(T.numBranches) //, nodes(T.nodes)
  {
    ids.resize(T.ids.size());
    for (size_t i(0); i < T.ids.size(); ++i) ids[i] = T.ids[i];
  }

  // Destructor
  virtual ~Tree() {}

  // ========================================================================
  // METHODS
  
  // read tree from file
  inline void read(string filename) { 
    fn = filename;
    // readTimes(fn,times,ttypes,extant,maxExtant); 
    readFromFile(fn);
  }
  void readFromFile(string fn);

  // get pointers
  inline double* ti() { return times.data(); }
  inline int* tt() { return ttypes.data(); }

  inline size_t size() const { return times.size(); }

  // set methods
  inline void set_nroot(int n) { nroot = n; }

  // get methods
  inline double time(size_t i) const { return times.at(i); }
  inline int ttype(size_t i) const { return ttypes.at(i); }
  inline int getMaxExtant() const { return maxExtant; }
  inline int getExtant() const { return extant; }

  inline size_t countTypes(int type) const {
    size_t c = 0;
    for (size_t i = 0; i < ttypes.size(); ++i) {
      if (ttype(i) == type) ++c;
    }
    return c;
  }

  inline int num_branches() const { return numBranches; }

  inline int max_id() const {
    int max = 0;
    for (size_t i = 0; i < ids.size(); ++i) {
      if (ids[i].size() > 0) {
        if (ids[i][0] > max) max = ids[i][0];
      }
    }
    return max;
  }

  // ========================================================================
  // OPERATORS

  const Tree& operator=(const Tree& t) {
    times = t.times;
    ttypes = t.ttypes;
    extant = t.extant;
    maxExtant = t.maxExtant;
    fn = t.fn;
    return *this;
  }

  // ========================================================================

  inline double maxTime() const { return times.back(); }
  void reverse();

  // ========================================================================

  void addRateShift(double t);

//  protected:
  // ========================================================================

  vector<double> times;     /* event times */
  vector<int> ttypes;       /* branching event type */
  int extant;               /* extant branches at t=0 */
  int maxExtant;            /* maximum number of extant branches */
  int nroot;                /* number of lineages at the root */
  string fn;                /* tree file name */
  int is_rev;               /* tree in forward or backwards format */
  // TreeNodeList nodes;    /* TreeNodes for extra information */
  int numBranches;
  vector< vector<int> > ids;
};

#endif
