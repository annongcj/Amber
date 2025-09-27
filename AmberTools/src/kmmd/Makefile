
AMBERSRC=../../..



INC   = -I$(AMBERSRC)/src/pmemd/src/cuda -I$(AMBERSRC)/src/pmemd/include -I$(AMBERTOOLS)/src/include 
DEFS  = -Duse_DPFP

#"-lpmemd_cuda_DPFP" "-lsander"

CUOBJS = kmmd_forcecalc.o hash.o pdb_io.o scanDB.o DB_vs_forcefield.o KMMDContext.o kmmd_externFortran.o KMMDParseJSON.o 

CLANG=-ccbin clang-3.8
CC=gcc
NVCC=nvcc
CFLAGS=-Wall
DBG=-D_DBG_UNITTESTS -g
NVCCFLAGS=$(_CLANG) $(DBG)

libkmmd.so: $(CUOBJS)
	$(NVCC) $(INC) $(DEFS) $(NVCCFLAGS) $(CUOBJS) -g -shared  -o libkmmd.so

testjson: KMMDParseJSON.cpp hash.o
	$(CC) $(INC) $(DEFS) -DDBG_JSON hash.c KMMDParseJSON.cpp -o testjson $(SHORTLINKFLAGS)

standalone: $(COBJS) $(CUOBJS) testwrap_forcecalc_cy.cu DB_snap_ff_eval.cpp
	$(NVCC) $(INC) $(DEFS) $(NVCCFLAGS) DB_snap_ff_eval.cpp testwrap_forcecalc_cy.cu $(CUOBJS) -o forcecalc $(LINKFLAGS)
	$(NVCC) $(INC) $(DEFS) -DTRAIN_STANDALONE $(NVCCFLAGS) DB_snap_ff_eval.cpp $(CUOBJS) -o dbsnap_ffeval $(LINKFLAGS)

%.o: %.c 
	$(NVCC) $(INC) $(DEFS) $(NVCCFLAGS) --compiler-options -fPIC -c -o $@ $<

%.o: %.cpp
	$(NVCC) $(INC) $(DEFS) $(NVCCFLAGS) --compiler-options -fPIC -c -o $@ $<

%.o: %.cu 
	$(NVCC) $(INC) $(DEFS) $(NVCCFLAGS) --compiler-options -fPIC -c -o $@ $<


.PHONY: clean
clean:
	rm -f *.o *~ \#*


