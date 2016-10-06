#ifndef SOLVER_HH
#define SOLVER_HH

class Solver {
protected:
  double m_dt;

public:
  Solver(double dt) : m_dt(dt) { }
  ~Solver() { }

  Solver(const Solver&) = delete;
  Solver &operator =(const Solver&) = delete;

  virtual void step(class SolarSystem &system) = 0;
};

class EulerSolver : public Solver {
public:
  EulerSolver(double dt) : Solver(dt) { }

  void step(class SolarSystem &system);
};

class VerletSolver : public Solver {
public:
  VerletSolver(double dt) : Solver(dt) { }

  void step(class SolarSystem &system);
};

#endif
