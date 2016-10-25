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
  double T  = 500;
  double dt = 0.01;

  // gravitational constant
  double G = 4 * M_PI * M_PI;

  // open file to write plot data to
  FILE *fp = fopen("3d.dat", "w");
  fprintf(fp, "v0\tt\tearth_x\tearth_y\tU\tK\n");

  // list of initial velocities (divided by 2π)
  std::vector<double> v0_list { 1.40, 1.41, 1.42, 1.43, 1.44 };

  for(double v0_k : v0_list) {
    double v0 = v0_k * 2 * M_PI;

    // create a solar system with the Earth and the Sun
    SolarSystem system(G);
    system.createFixedBody(vec3(0, 0, 0), 1);               // Sun
    system.createBody(vec3(1, 0, 0), vec3(0, v0, 0), 3e-6); // Earth
    CelestialBody &sun   = system.body(0);
    CelestialBody &earth = system.body(1);

    // run simulation
    double t;
    int i;
    VerletSolver solver;
    system.update();
    for(t = 0, i = 0; t <= T; i++) {
      // for every hundredth point, write solar system state to file
      if(i % 100 == 0) {
        fprintf(fp, "%.2f\t%.4E\t%.4E\t%.4E\t%.4E\t%.4E\n",
          v0_k, t,
          earth.position.x(), earth.position.y(),
          system.potentialEnergy(),
          system.kineticEnergy()
        );
      }

      solver.step(system, dt);
      t += dt;
    }
  }
  fclose(fp);
}
