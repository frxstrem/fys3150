CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

.PHONY: all

all: b b.test d

clean:
	rm -rf b d *.o *.dat

b: b.cc comp_eig.hh
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}

b.test: b.test.cc b.cc comp_eig.hh
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}

d: d.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${LIBS}
