#include "solver.hh"
#include "solarsystem.h"
#include "vec3.h"

void EulerSolver::step(class SolarSystem &system) {
  system.calculateForcesAndEnergy();

  for(CelestialBody &body : system.bodies()) {
    body.velocity += body.force / body.mass * m_dt;
    body.position += body.velocity * m_dt;
  }
}

void VerletSolver::step(class SolarSystem &system) {
  system.calculateForcesAndEnergy();

  for(CelestialBody &body : system.bodies()) {
    vec3 acc = body.force / body.mass;

    body.position += body.velocity * m_dt + 0.5 * acc * m_dt * m_dt;
    body.velocity += acc * m_dt / 2;
  }

  system.calculateForcesAndEnergy();

  for(CelestialBody &body : system.bodies()) {
    vec3 acc = body.force / body.mass;
    body.velocity += acc * m_dt / 2;
  }
}
