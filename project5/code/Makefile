CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

MPI_CXX=mpic++

.PHONY: all

all: 5a 5a-equilibrium 5c 5d 5e

5a: 5a.o transactions.o histograms.o
	${MPI_CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

5a-equilibrium: 5a-equilibrium.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

5c: 5c.o transactions.o histograms.o
	${MPI_CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

5d: 5d.o transactions.o histograms.o
	${MPI_CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

5e: 5e.o transactions.o histograms.o
	${MPI_CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

5a.o: 5a.cc transactions.hh histograms.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<

5a-equilibrium.o: 5a-equilibrium.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<

5c.o: 5c.cc transactions.hh histograms.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<

5d.o: 5d.cc transactions.hh histograms.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<

5e.o: 5e.cc transactions.hh histograms.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<

transactions.o: transactions.cc transactions.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<

histograms.o: histograms.cc histograms.hh
	${MPI_CXX} ${CXXFLAGS} -c -o $@ $<
