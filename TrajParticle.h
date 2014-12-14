#ifndef __TRAJPARTICLE_H__
#define __TRAJPARTICLE_H__

#include <stdlib.h>
#include <unistd.h>

#include "MCTraj.h"
#include "Trajectory.h"
#include "TrajParticle.h"
#include "ParticleFilter.h"
#include "../Tree.h"

#include <gsl/gsl_rng.h>

namespace MCTraj {
  class TrajParticle : public Trajectory, public Particle {
    public:
      TrajParticle() 
        : Trajectory(), Particle(), init_time(0.0) {}

      TrajParticle(string name, double weight = 1.0) 
        : Trajectory(), Particle(name,weight), init_time(0.0)
      {}

      TrajParticle(string name, double weight, const Trajectory& traj)
        : Trajectory(traj), Particle(name,weight), init_time(0.0)
      {}

      TrajParticle(const TrajParticle& tp) 
        : Trajectory(tp), Particle(tp), init_time(tp.init_time), ivals(tp.ivals)
      {}

      virtual ~TrajParticle() {}

      const TrajParticle& copy(const TrajParticle& tp) {
        Trajectory::copyState(tp);
        Particle::operator=(tp);
        init_time = tp.time;
        ivals = tp.ivals;
        return *this;
      }

      const TrajParticle& operator=(const TrajParticle& tp) { 
        Trajectory::operator=(tp);
        Particle::operator=(tp);
        init_time = tp.init_time;
        ivals = tp.ivals;
        return *this;
      }

      void setLast() { ivals.push_back(transitions.size()-1); }
      // void setLast(int l) { last = l; }

      int getLast() const { return (ivals.size() > 0) ? ivals.back() : -1; }
      double getInit() const { return init_time; }

      int lastType(int type) const;

      void resetWeight() {
        Particle::resetWeight();
        prob = 1.0;
      }

      void print_from_last(ostream& out, string prefix = "", const void* pars = NULL, int offset = 0) const;
      void print_from_first(ostream& out, string prefix = "", const void* pars = NULL) const;

    protected:
      double init_time;
      // int last;
      vector<int> ivals;
  };
}

#endif

