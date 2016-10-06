#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <memory>
#include <stdlib.h>
#include "solarsystem.h"
#include "solver.hh"
using namespace std;

struct SimulationParameters {
  // solar system object
  shared_ptr<SolarSystem> system;

  // number of time steps
  int M;
  // length of time steps
  double dt;
};

static SimulationParameters read_params_file(const char *filename);

/**
 * Read simulation parameters from file (uses SolarSystem.dat if none is specified),
 * then run simulation.
**/
int main(int argc, char **argv)
{
  const char *filename;
  if(argc > 1)
    filename = argv[1];
  else
    filename = "SolarSystem.dat";

  // read simulation parameters from file
  SimulationParameters params = read_params_file(filename);
  SolarSystem &system = *params.system;
  int M = params.M;
  double dt = params.dt;

  // print list of celestial bodies
  cout << "Celestial bodies:" << endl;
  for(CelestialBody &body : system.bodies())
    cout << "  mass=" << body.mass << ", pos=" << body.position << ", vel=" << body.velocity << endl;

  // integrate
  EulerSolver solver(dt);
  for(int m = 0; m < M; m++) {
    solver.step(system);
    system.writeToFile("positions.xyz");
  }

  cout << "DONE!" << endl;

  // print list of celestial bodies
  cout << "Celestial bodies:" << endl;
  for(CelestialBody &body : system.bodies())
    cout << "  mass=" << body.mass << ", pos=" << body.position << ", vel=" << body.velocity << endl;
}

static SimulationParameters read_params_file(const char *filename) {
  SimulationParameters params;
  ifstream in(filename);

  // IMPORTANT! this code *assumes* that the file is in the right format

  // first line: comment line (ignore)
  string line; getline(in, line);

  // second line: total simulation time
  // third line: simulation time steps
  double T, dt;
  in >> T >> dt;

  // calculate number of steps M
  params.M = T / dt;
  params.dt = dt;

  // fourth line: gravitational constant (divided by 4 * π²)
  double G;
  in >> G;
  G *= 4 * M_PI * M_PI;

  // fifth line: number of planets N
  int N;
  in >> N;

  // create solar system
  params.system = make_shared<SolarSystem>(G);

  // next N lines: (m x y z vx vy vz) initial conditions for each body
  for(int i = 0; i < N; i++) {
    double m, x, y, z, vx, vy, vz;
    in >> m >> x >> y >> z >> vx >> vy >> vz;

    params.system->createCelestialBody( vec3(x, y, z), vec3(vx, vy, vz), m );
  }

  return params;
}
