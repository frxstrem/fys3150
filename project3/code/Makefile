CXX=g++
CXXFLAGS=--std=c++14 -O0
LDFLAGS=
LIBS=-lm

.PHONY: all clean

all: program

clean:
	@rm -rf program *.o *.xyz

program: main.o celestialbody.o solver.o solarsystem.o vec3.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

main.o: main.cpp solarsystem.h solver.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

solarsystem.o: solarsystem.cpp celestialbody.h solarsystem.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

celestialbody.o: celestialbody.cpp celestialbody.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

solver.o: solver.cc solver.hh vec3.h
	${CXX} ${CXXFLAGS} -c -o $@ $<

vec3.o: vec3.cpp vec3.h
	${CXX} ${CXXFLAGS} -c -o $@ $<
