# Makefile for Gaussian accelerated molecular dynamics tests

all: test

pme:
	@echo "Testing GaMD with PME"
	@-cd dhfr_pme && ./Run.pme.gamd1 $(OPT)
	@-cd dhfr_pme && ./Run.pme.gamd2 $(OPT)
	@-cd dhfr_pme && ./Run.pme.gamd3 $(OPT)

gb:
	@echo "Testing GaMD with GB"
	@-cd rna_gb && ./Run.gb.gamd1 $(OPT)
	@-cd rna_gb && ./Run.gb.gamd2 $(OPT)
	@-cd rna_gb && ./Run.gb.gamd3 $(OPT)

test: pme gb

pmemd: 
	-$(MAKE) test OPT=pmemd
