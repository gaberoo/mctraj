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

class TreeBranch {
  public:
    TreeBranch() : id(-1), type(-1), time(0.0), child1(-1), child2(-1) {}
    TreeBranch(int id, int type, int time, int c1 = -1, int c2 = -1)
      : id(id), type(type), time(time), child1(c1), child2(c2)
    {}
    TreeBranch(const TreeBranch& b)
      : id(b.id), type(b.type), time(b.time), 
        child1(b.child1), child2(b.child2)
    {}
    virtual ~TreeBranch() {}

    int id;
    int type;
    double time;
    int child1;
    int child2;
};

class Tree {
public:
  /* CONSTRUCTORS-DESTRUCTORS */

  // 1. Default
  Tree() 
    : extant(0), maxExtant(0), nroot(0), fn(""), 
      is_rev(0), numBranches(0) 
  {}

  // 2. Read tree from file
  Tree(const char* fn, int nroot = 0) 
    : nroot(nroot), is_rev(0)
  { read(string(fn)); }

  Tree(string fn, int nroot = 0) 
    : nroot(nroot), is_rev(0)
  { read(fn); }

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
  inline int read(string filename) { 
    fn = filename;
    // readTimes(fn,times,ttypes,extant,maxExtant); 
    return readFromFile(fn);
  }
  int readFromFile(string fn);

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

  size_t countTypes(int type) const;
  inline int num_branches() const { return numBranches; }
  int max_id() const;

  // ========================================================================
  // OPERATORS

  inline const Tree& operator=(const Tree& t) {
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
  void makeBranches();

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
  vector<TreeBranch> branches;
};

#endif
