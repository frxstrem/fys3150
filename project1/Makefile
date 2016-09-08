CXX=g++
CXXFLAGS=--std=c++17
LDFLAGS=-lm
LIBS=

.PHONY: all list clean

all: main pdf

list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

clean:
	@git clean -fxd .

pdf: project1.tex
	latexmk -pdf -interaction=nonstopmode -pdflatex="lualatex -file-line-error --shell-escape" project1.tex

main: main.o
	${CXX} ${LDFLAGS} -o $@ $^ ${LIBS}

main.o: main.cc mat/Matrix.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<
