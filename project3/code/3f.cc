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
  double T  = 25;
  double dt = 0.001;

  // gravitational constant
  double G = 4 * M_PI * M_PI;

  // open file to write plot data to
  FILE *fp = fopen("3f.dat", "w");
  fprintf(fp, "body\tt\tx\ty\tz\n");

  // NASA's velocities are expressed as AU/day, while the simulation uses AU/year
  // `s' is the conversion factor
  double s = 365.24;

  // create a solar system with all the planets + the Sun and the Moon
  SolarSystem system(G);
  system.createBody( // Sun
    vec3( 3.547262250463157E-03,  3.471009732655305E-03, -1.594233907062380E-04),
    vec3(-2.073074446530653E-06,  6.815546336340393E-06,  4.271186269439737E-08) * s,
    1);
  system.createBody( // Mercury
    vec3(-3.840022217498515E-01, -1.571218382135171E-01,  2.227270778117466E-029),
    vec3( 4.951859715755229E-03, -2.476481319335340E-02, -2.478579232708704E-03) * s,
    3.3e23 / 2e30
  );
  system.createBody( // Venus
    vec3( 3.669723300842595E-01, -6.267508777559652E-01, -2.977275548926195E-02),
    vec3( 1.738419265724836E-02,  1.004173699067179E-02, -8.656583757242518E-04) * s,
    4.9e24 / 2e30
  );
  system.createBody( // Earth
    vec3( 8.580832100932452E-01,  5.124323971262588E-01, -1.814930137100218E-04),
    vec3(-9.081258370972716E-03,  1.472861977023750E-02, -1.023214727329336E-06) * s,
    6e24 / 2e30
  );
  system.createBody( // Moon
    vec3( 8.562465483701344E-01,  5.142498654067389E-01, -2.839390710657199E-04),
    vec3(-9.515827954093179E-03,  1.433997381647825E-02,  4.466403021782183E-05) * s,
    7.35e22 / 2e30
  );
  system.createBody( // Mars
    vec3( 1.230144669794901E+00, -6.305833129893422E-01, -4.354973947040395E-02),
    vec3( 6.956734020408076E-03,  1.363624768289644E-02,  1.148717471176532E-04) * s,
    6.6e23 / 2e30
  );
  system.createBody( // Jupiter
    vec3(-5.423289894366403E+00, -5.180084449201873E-01,  1.234383246843473E-01),
    vec3( 6.294052333429065E-04, -7.155995867542226E-03,  1.560106895550821E-05) * s,
    1.9e27 / 2e30
  );
  system.createBody( // Saturn
    vec3(-2.220861899049674E+00, -9.786202215943753E+00,  2.585434814126087E-01),
    vec3( 5.134758705420029E-03, -1.251502185262647E-03, -1.829320865904290E-04) * s,
    5.5e26 / 2e30
  );
  system.createBody( // Uranus
    vec3( 1.844952322607655E+01,  7.592520781423359E+00, -2.108185350546166E-01),
    vec3(-1.525630336541416E-03,  3.453808892579247E-03,  3.246304225728262E-05) * s,
    8.8e25 / 2e30
  );
  system.createBody( // Neptune
    vec3( 2.827009173035302E+01, -9.895464694829448E+00, -4.477347883067835E-01),
    vec3( 1.016213511636345E-03,  2.981128557772945E-03, -8.517801554511235E-05) * s,
    1.03e26 / 2e30
  );

  // give the sun velocity such that total momentum is zero
  CelestialBody &sun = system.body(0);
  vec3 momentum;
  for(CelestialBody &body : system.bodies())
    momentum += body.mass * body.velocity;
  sun.velocity -= momentum / sun.mass;

  // run simulation
  double t;
  int i;
  VerletSolver solver;
  system.update();
  for(t = 0, i = 0; t <= T; i++) {
    // for every 25th point, write solar system state to file
    if(i % 25 == 0) {
      for(int j = 0; j < system.numOfBodies(); j++) {
        CelestialBody &body = system.body(j);

        fprintf(fp, "%d\t%.4E\t%.4E\t%.4E\t%.4E\n",
          j, t,
          body.position.x(), body.position.y(), body.position.z()
        );
      }
    }

    solver.step(system, dt);
    t += dt;
  }
  fclose(fp);
}
