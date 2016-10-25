CXX=g++
CXXFLAGS=--std=c++14
LDFLAGS=
LIBS=-lm

.PHONY: all clean

all: 3b 3d 3e 3f 3g

clean:
	@rm -rf *.o

3b: 3b.o CelestialBody.o Solver.o SolarSystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

3b.o: 3b.cc SolarSystem.hh Solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

3d: 3d.o CelestialBody.o Solver.o SolarSystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

3d.o: 3d.cc SolarSystem.hh Solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

3e: 3e.o CelestialBody.o Solver.o SolarSystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

3e.o: 3e.cc SolarSystem.hh Solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

3f: 3f.o CelestialBody.o Solver.o SolarSystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

3f.o: 3f.cc SolarSystem.hh Solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

3g: 3g.o CelestialBody.o Solver.o SolarSystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

3g.o: 3g.cc SolarSystem.hh Solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

SolarSystem.o: SolarSystem.cc CelestialBody.hh SolarSystem.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

CelestialBody.o: CelestialBody.cc CelestialBody.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

Solver.o: Solver.cc Solver.hh vec3.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

vec3.o: vec3.cc vec3.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<
