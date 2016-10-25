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
  double T  = 100;
  double dt = 0.000001;

  // gravitational constant
  double G = 4 * M_PI * M_PI;

  // open file to write plot data to
  FILE *fp = fopen("3g.dat", "w");
  fprintf(fp, "t\tx\ty\tangle\n");

  // NASA's velocities are expressed as AU/day, while the simulation uses AU/year
  // `s' is the conversion factor
  double s = 365.24;

  // create a solar system with Mercury and the Sun
  SolarSystem system(G);
  system.createFixedBody(vec3(0, 0, 0), 1);
  system.createBody(vec3(0.3075, 0, 0), vec3(0, 12.44, 0), 3.3e23 / 2e30);

  CelestialBody &sun     = system.body(0);
  CelestialBody &mercury = system.body(1);

  // run simulation
  double t;
  int i;
  VerletSolver solver;
  system.update();

  for(t = 0, i = 0; t <= T; i++) {
    if(i % int(0.1 / dt) == 0)
      printf("t = %5.1f\n", t);

    if(mercury.velocity.length() >= 12.44 - 1e-10) {
      fprintf(fp, "%.4E\t%.4E\t%.4E\t%.4E\n",
        t,
        mercury.position.x(), mercury.position.y(),
        atan2(mercury.position.y(), mercury.position.x()) * (180 / M_PI)
      );
    }

    solver.step(system, dt);
    t += dt;
  }
  fclose(fp);
}
