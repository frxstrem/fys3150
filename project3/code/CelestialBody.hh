#ifndef CELESTIAL_BODY_HH
#define CELESTIAL_BODY_HH

#include "vec3.hh"

class CelestialBody {
public:
  vec3 position;
  vec3 velocity;
  vec3 acceleration;
  double mass;
  const bool fixed;

  CelestialBody(vec3 position, vec3 velocity, double mass)
    : position(position), velocity(velocity), acceleration(), mass(mass), fixed(false)
  { }
  CelestialBody(vec3 position, vec3 velocity, double mass, bool fixed)
    : position(position), velocity(velocity), acceleration(), mass(mass), fixed(fixed)
  { }

  CelestialBody(const CelestialBody &body)
    : position(body.position), velocity(body.velocity), acceleration(body.acceleration), mass(body.mass), fixed(body.fixed)
  { }

  // is this celestial body fixed to a specific position?
  bool isFixed() const { return fixed; }
};

#endif
