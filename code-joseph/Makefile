CXX=g++
CXXFLAGS=--std=c++14 -Ofast
LDFLAGS=
LIBS=-lm -larmadillo

.PHONY: list all

list:
	@echo "Commands:"
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

all: jacobi jacobi_step_data

clean:
	rm b *.o

jacobi: jacobi.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

jacobi_step_data: jacobi_step_data.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<
