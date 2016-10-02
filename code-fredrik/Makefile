CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

.PHONY: all

all: b d

clean:
	rm -rf b d *.o *.dat

b: b.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}
d: d.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<
