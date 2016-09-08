CXX=g++
CXXFLAGS=--std=c++17
LDFLAGS=-lm
LIBS=

.PHONY: build clean

build: main

clean:
	@rm -rf main *.o

main: main.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

main.o: main.cc Matrix.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<
