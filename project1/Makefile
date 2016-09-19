CXX=g++
CXXFLAGS=--std=c++17
LDFLAGS=-lm
LIBS=

.PHONY: list all run clean

list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

all: run pdf

clean:
	@git clean -fxd .

pdf: project1.tex
	latexmk -pdf -interaction=nonstopmode -pdflatex="lualatex -file-line-error --shell-escape" project1.tex

run: oppg_b oppg_d oppg_e
	./oppg_b
	./oppg_d
	./oppg_e

oppg_b: oppg_b.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

oppg_d: oppg_d.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

oppg_e: oppg_e.o lib.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

%.o: %.cc
	${CXX} ${CXXFLAGS} -c -o $@ $<
