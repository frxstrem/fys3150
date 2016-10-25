#include <cstdio>
#include "vec3.hh"
#include "SolarSystem.hh"
#include "CelestialBody.hh"
#include "Solver.hh"

void EulerSolver::step(class SolarSystem &system, double dt) const {
  system.update();

  for(CelestialBody &body : system.bodies()) {
    if(body.isFixed()) continue;

    body.position += body.velocity * dt;
    body.velocity += body.acceleration * dt;
  }
}

void VerletSolver::step(class SolarSystem &system, double dt) const {
  system.update();

  for(CelestialBody &body : system.bodies()) {
    if(body.isFixed()) continue;

    body.position += body.velocity * dt + 0.5 * body.acceleration * dt * dt;
    body.velocity += 0.5 * body.acceleration * dt;
  }

  system.update();

  for(CelestialBody &body : system.bodies()) {
    if(body.isFixed()) continue;

    body.velocity += 0.5 * body.acceleration * dt;
  }
}
