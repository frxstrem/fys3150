CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

.PHONY: list all

list:
	@echo "Commands:"
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

all: b

clean:
	rm b *.o

b: b.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<
