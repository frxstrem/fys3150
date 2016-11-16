CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

MPI_CXX=mpic++

.PHONY: all

all: 4b 4c 4d 4e

4b: 4b.cc IsingState.hh metropolis.hh misc.hh
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}

4c: 4c.cc IsingState.hh metropolis.hh misc.hh
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}

4d: 4d.cc IsingState.hh metropolis.hh misc.hh
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}

4e: 4e.cc IsingState.hh metropolis.hh misc.hh
	${MPI_CXX} ${CXXFLAGS} ${MPI_CFLAGS} ${LDFLAGS} -o $@ $< ${LIBS} ${MPI_LIBS}
