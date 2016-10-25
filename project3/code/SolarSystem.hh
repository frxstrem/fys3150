#ifndef SOLAR_SYSTEM_HH
#define SOLAR_SYSTEM_HH

#include <vector>
#include "vec3.hh"
#include "CelestialBody.hh"

class SolarSystem {
protected:
  double mG; // gravitational constant

private:
  std::vector<CelestialBody> mBodies;

  double mKineticEnergy, mPotentialEnergy;
  vec3 mTotalAngularMomentum;

public:
  SolarSystem(double G)
    : mG(G), mBodies(), mKineticEnergy(),
      mPotentialEnergy(), mTotalAngularMomentum()
  { }

  // create celestial body
  void createBody(vec3 position, vec3 velocity, double mass);
  void createFixedBody(vec3 position, double mass);

  // get number of celestial bodies
  int numOfBodies() const;

  // get a specific celestial body by index
  CelestialBody &body(int index);

  // get iterator of all celestial bodies
  std::vector<CelestialBody> &bodies();
  const std::vector<CelestialBody> &bodies() const;

  // calculate gravitational force acting on first body
  virtual vec3 calcForce(const CelestialBody &a, const CelestialBody &b) const;

  // recalculate forces, energies and total angular momentum
  void update();

  // get energies
  double kineticEnergy() const { return mKineticEnergy; }
  double potentialEnergy() const { return mPotentialEnergy; }
  double totalEnergy() const { return kineticEnergy() + potentialEnergy(); }

  // get angular momentum
  vec3 angularMomentum() const { return mTotalAngularMomentum; }
};

class RelativisticSolarSystem : public SolarSystem {
public:
  RelativisticSolarSystem(double G) : SolarSystem(G) { }

  virtual vec3 calcForce(const CelestialBody &a, const CelestialBody &b) const;
};

#endif
