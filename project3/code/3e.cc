#include <cmath>

#include "vec3.hh"
#include "SolarSystem.hh"
#include "CelestialBody.hh"
#include "Solver.hh"

/**
 * Units are expressed in astronomical units (distance), years (time) and solar masses (mass).
**/

int main(int argc, char **argv)
{
  // set simulation duration and time step
  double T  = 5;
  double dt = 0.001;

  // gravitational constant
  double G = 4 * M_PI * M_PI;

  // open file to write plot data to
  FILE *fp = fopen("3e.dat", "w");
  fprintf(fp, "t\tearth_x\tearth_y\tearth_vx\tearth_vy\tsun_x\tsun_y\tjupiter_x\tjupiter_y\tK\tU\n");

  // create a solar system with the Earth and the Sun
  SolarSystem system(G);
  system.createFixedBody(vec3(0, 0, 0), 1);                                     // Sun
  system.createBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6);                 // Earth
  system.createBody(vec3(5.20, 0, 0), vec3(0, 2 * M_PI / sqrt(5.20), 0), 1e-3); // Jupiter

  CelestialBody &sun     = system.body(0);
  CelestialBody &earth   = system.body(1);
  CelestialBody &jupiter = system.body(2);

  // run simulation
  double t;
  int i;
  VerletSolver solver;
  system.update();
  for(t = 0, i = 0; t <= T; i++) {
    // for every tenth point, write solar system state to file
    if(i % 10 == 0) {
      fprintf(fp, "%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\n",
        t,
        earth.position.x(), earth.position.y(),
        earth.velocity.x(), earth.velocity.y(),
        sun.position.x(), sun.position.y(),
        jupiter.position.x(), jupiter.position.y(),
        system.kineticEnergy(),
        system.potentialEnergy()
      );
    }

    solver.step(system, dt);
    t += dt;
  }
  fclose(fp);
}
