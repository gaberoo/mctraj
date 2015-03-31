#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <iostream>
#include <cmath>
using namespace std;

namespace MCTraj {
  class Particle {
    public:
      Particle() {}

      Particle(string n, double w) 
        : name(n), weight(w), cum_weight(0.0), id(0), parent(0), msg("")
      {}

      Particle(const Particle& p) 
        : name(p.name), weight(p.weight), cum_weight(p.cum_weight), 
          id(p.id), parent(p.parent), msg(p.msg)
      {}

      virtual ~Particle() {}

      inline double getWeight() const { return weight; }
      inline string getName() const { return name; }
      inline int getId() const { return id; }
      inline int getParent() const { return parent; }

      inline void setWeight(double w) { weight = w; }
      inline void updateWeight(double w) { weight *= w; }

      inline void setName(string n) { name = n; }
      inline void setId(int i) { id = i; }
      inline void setParent(int p) { parent = p; }

      const Particle& operator=(const Particle& p) { 
        name = p.name;
        weight = p.weight; 
        cum_weight = p.cum_weight;
        id = p.id;
        parent = p.parent;
        msg = p.msg;
        return *this;
      }

      virtual void resetWeight() {
        cum_weight += log(weight);
        weight = 1.0;
      }

      inline double cw() const { return cum_weight + log(weight); }

    protected:
      string name;
      double weight;
      double cum_weight;
      int id;
      int parent;
      string msg;
  };
}

#endif

