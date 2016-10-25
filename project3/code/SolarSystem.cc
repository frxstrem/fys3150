#include <vector>

#include "vec3.hh"
#include "CelestialBody.hh"
#include "SolarSystem.hh"

void SolarSystem::createBody(vec3 position, vec3 velocity, double mass) {
  CelestialBody body(position, velocity, mass);
  mBodies.push_back(body);
}

void SolarSystem::createFixedBody(vec3 position, double mass) {
  CelestialBody body(position, vec3(), mass, true);
  mBodies.push_back(body);
}

int SolarSystem::numOfBodies() const {
  return mBodies.size();
}

CelestialBody &SolarSystem::body(int index) {
  return mBodies[index];
}

std::vector<CelestialBody> &SolarSystem::bodies() {
  return mBodies;
}

const std::vector<CelestialBody> &SolarSystem::bodies() const {
  return mBodies;
}

vec3 SolarSystem::calcForce(const CelestialBody &a, const CelestialBody &b) const {
  vec3 deltaPosition = a.position - b.position;
  double dr = deltaPosition.length();
  return -(mG * a.mass * b.mass / (dr * dr)) * (deltaPosition / dr);
}

vec3 RelativisticSolarSystem::calcForce(const CelestialBody &a, const CelestialBody &b) const {
  vec3 deltaPosition = a.position - b.position;
  vec3 deltaVelocity = a.velocity - b.velocity;
  double dr = deltaPosition.length();

  // speed of light (in AU/yr)
  static const double c = 63198;

  // calculate angular momentum per unit mass
  double l = deltaPosition.cross(deltaVelocity).length();

  // calculate relativistic gravitational force
  double k = l / (dr * c);
  return -(mG * a.mass * b.mass / (dr * dr)) * (1 + 3 * k * k) * (deltaPosition / dr);
}

void SolarSystem::update() {
  mKineticEnergy = 0;
  mPotentialEnergy = 0;
  mTotalAngularMomentum.zeros();

  for(CelestialBody &body : mBodies) {
    body.acceleration.zeros();

    if(body.isFixed())
      body.velocity.zeros();
  }

  for(int i = 0; i < mBodies.size(); i++) {
    CelestialBody &A = mBodies[i];   // first body
    for(int j = i + 1; j < mBodies.size(); j++) {
      CelestialBody &B = mBodies[j]; // second body

      // calculate force and potential energy
      vec3 F = calcForce(A, B);
      A.acceleration += F / A.mass;
      B.acceleration -= F / B.mass;

      mPotentialEnergy -= mG * A.mass * B.mass / (A.position - B.position).length();
    }

    // calculate kinetic energy and angular momentum
    mKineticEnergy += 0.5 * A.mass * A.velocity.lengthSquared();
    mTotalAngularMomentum += A.position.cross(A.mass * A.velocity);
  }
}
