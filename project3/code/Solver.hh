#ifndef SOLVER_HH
#define SOLVER_HH

#include "vec3.hh"
#include "SolarSystem.hh"

class Solver {
public:
  Solver() { }
  ~Solver() { }

  Solver(const Solver&) = delete;
  Solver &operator =(const Solver&) = delete;

  virtual void step(class SolarSystem &system, double dt) const = 0;
};

class EulerSolver: public Solver {
public:
  virtual void step(class SolarSystem &system, double dt) const;
};

class VerletSolver: public Solver {
public:
  virtual void step(class SolarSystem &system, double dt) const;
};

#endif
