.PHONY: list all run clean

list:
	@echo "Commands:"
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' | xargs

all: run pdf

clean:
	@git clean -fxd .

pdf: project2.tex
	latexmk -pdf -interaction=nonstopmode -pdflatex="lualatex -file-line-error --shell-escape" project2.tex

run: programs

programs:
	make -C code all
