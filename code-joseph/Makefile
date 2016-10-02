CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

.PHONY: all

all: jacobi jacobi_step_data jacobi_interact

clean:
	rm -rf jacobi jacobi_step_data jacobi_interact *.o *.dat

jacobi: jacobi.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

jacobi_step_data: jacobi_step_data.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

jacobi_interact: jacobi_interact.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<
