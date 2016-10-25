#include <cmath>

#include "vec3.hh"
#include "SolarSystem.hh"
#include "CelestialBody.hh"
#include "Solver.hh"

/**
 * Units are expressed in astronomical units (distance), years (time) and solar masses (mass).
**/

void run_simulation(const char *filename, const Solver &solver, double T, double dt);

int main(int argc, char **argv)
{
  // set simulation duration and time step
  double T  = 5;
  double dt = 0.001;

  // run simulation with Euler's forward algorithm
  run_simulation("3b_euler.dat", EulerSolver(), T, dt);

  // run simulation with Verlet algorithm
  run_simulation("3b_verlet.dat", VerletSolver(), T, dt);
}

void run_simulation(const char *filename, const Solver &solver, double T, double dt) {
  double G = 4 * M_PI * M_PI; // gravitational constant

  // create a solar system with the Earth and the Sun
  SolarSystem system(G);
  system.createFixedBody(vec3(0, 0, 0), 1);                     // Sun
  system.createBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6); // Earth
  CelestialBody &sun   = system.body(0);
  CelestialBody &earth = system.body(1);

  // open file to write plot data to
  // columns:
  //   time    earth x,y position    sun x,y position    potential and kinetic energy    total angular momentum
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "t\tearth_x\tearth_y\tsun_x\tsun_y\tU\tK\tLx\tLy\tLz\n");

  // run simulation
  double t;
  int i;
  system.update();
  for(t = 0, i = 0; t <= T; i++) {
    // for every tenth point, write solar system state to file
    if(i % 10 == 0) {
      fprintf(fp, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\n",
        t,
        earth.position.x(), earth.position.y(),
        sun.position.x(), sun.position.y(),
        system.potentialEnergy(),
        system.kineticEnergy(),
        system.angularMomentum().x(),
        system.angularMomentum().y(),
        system.angularMomentum().z()
      );
    }

    solver.step(system, dt);
    t += dt;
  }
  fclose(fp);
}
