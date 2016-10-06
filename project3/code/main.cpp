#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "solver.hh"
using namespace std;

int main(int numArguments, char **arguments)
{
  int numTimesteps = 1000;
  if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

  SolarSystem solarSystem(4 * M_PI * M_PI); // pass gravitational constant as a paremeter

  CelestialBody &sun = solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );

  // We don't need to store the reference, but just call the function without a left hand side
  solarSystem.createCelestialBody( vec3(1.5, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6 );

  // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
  vector<CelestialBody> &bodies = solarSystem.bodies();

  for(int i = 0; i<bodies.size(); i++) {
    CelestialBody &body = bodies[i]; // Reference to this body
    cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
  }

  double dt = 0.001;
  EulerSolver integrator(dt);
  for(int timestep=0; timestep<numTimesteps; timestep++) {
    integrator.step(solarSystem);
    solarSystem.writeToFile("positions.xyz");
  }

  cout << "I just created my first solar system that has " << solarSystem.bodies().size() << " objects." << endl;
  return 0;
}
