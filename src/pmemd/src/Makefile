#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999, 2004, 2008 **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************

# Makefile for Amber 24 PMEMD
SHELL=/bin/sh
AMBER=Amber24
CONFIG_FILE=../../config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA

# Platform-specific info should be found in ../../config.h
include $(CONFIG_FILE)

# PLUMED
ifdef PLUMED_INCLUDE_FILE
include $(PLUMED_INCLUDE_FILE)
endif

.NOTPARALLEL: cuda_SPFP cuda_DPFP pmemd.cuda_SPFP$(SFX) pmemd.cuda_DPFP$(SFX)

OBJS=   random$(OSFX) resamplekin$(OSFX) gbl_constants$(OSFX) \
        gbl_datatypes$(OSFX) state_info$(OSFX) file_io_dat$(OSFX) \
        boost-gamma-distribution$(OSFX) external_dat$(OSFX) \
        mdin_ctrl_dat$(OSFX) mdin_emil_dat$(OSFX) mdin_ewald_dat$(OSFX) mdin_debugf_dat$(OSFX) prmtop_dat$(OSFX) \
        inpcrd_dat$(OSFX) dynamics_dat$(OSFX) emil$(OSFX) img$(OSFX) nbips$(OSFX) offload_allocation$(OSFX) \
        parallel_dat$(OSFX) parallel$(OSFX) gb_parallel$(OSFX) \
        pme_direct$(OSFX) pme_recip_dat$(OSFX) pme_slab_recip$(OSFX) pme_blk_recip$(OSFX) \
        pme_slab_fft$(OSFX) pme_blk_fft$(OSFX) pme_fft_dat$(OSFX) fft1d$(OSFX) \
        bspline$(OSFX) nebread$(OSFX) neb$(OSFX) pme_force$(OSFX) pbc$(OSFX) nb_pairlist$(OSFX) gb_ene_hybrid$(OSFX) \
        nb_exclusions$(OSFX) cit$(OSFX) dynamics$(OSFX) bonds$(OSFX) angles$(OSFX) dihedrals$(OSFX) \
        extra_pnts_nb14$(OSFX) runmd$(OSFX) loadbal$(OSFX) shake$(OSFX) prfs$(OSFX) mol_list$(OSFX) \
        runmin$(OSFX) constraints$(OSFX) axis_optimize$(OSFX) gb_ene$(OSFX) veclib$(OSFX) phmd$(OSFX) gb_force$(OSFX) \
        timers$(OSFX) pmemd_lib$(OSFX) runfiles$(OSFX) file_io$(OSFX) \
        AmberNetcdf$(OSFX) bintraj$(OSFX) binrestart$(OSFX) pmemd_clib$(OSFX) \
        pmemd$(OSFX) degcnt$(OSFX) erfcfun$(OSFX) nmr_calls$(OSFX) nmr_lib$(OSFX) \
        get_cmdline$(OSFX) master_setup$(OSFX) pme_alltasks_setup$(OSFX) pme_setup$(OSFX) \
        ene_frc_splines$(OSFX) gb_alltasks_setup$(OSFX) nextprmtop_section$(OSFX) \
        angles_ub$(OSFX) dihedrals_imp$(OSFX) cmap$(OSFX) charmm$(OSFX) charmm_gold$(OSFX) \
        findmask$(OSFX) remd$(OSFX) multipmemd$(OSFX) remd_exchg$(OSFX) amd$(OSFX) gamd$(OSFX) ti$(OSFX) gbsa$(OSFX) \
        barostats$(OSFX) scaledMD$(OSFX) constantph$(OSFX) constante$(OSFX) energy_records$(OSFX) constantph_dat$(OSFX) \
        constante_dat$(OSFX) external$(OSFX) relaxmd$(OSFX) sgld$(OSFX) emap$(OSFX) get_efield_energy$(OSFX) \
        gti$(OSFX) reaf$(OSFX) \
        processor$(OSFX) parallel_processor$(OSFX) pme_fft_midpoint$(OSFX) pme_recip_midpoint$(OSFX) \
        dihedrals_midpoint$(OSFX) angles_midpoint$(OSFX) bonds_midpoint$(OSFX) dbg_arrays$(OSFX) mcres$(OSFX) \
        asm_CV$(OSFX) asm_multiCV$(OSFX) asm_splines_utilities$(OSFX) asm_utilities$(OSFX) asm_CV_utilities$(OSFX) asm$(OSFX) \
        ensure_alloc$(OSFX) reservoir${OSFX} sams${OSFX} ramd${OSFX} runreweight${OSFX} hybridsolvent_remd$(OSFX) \
		md_scheme$(OSFX)

NFE_OBJECTS = \
        nfe_lib.o nfe_setup.o nfe_colvar.o nfe_smd.o nfe_abmd.o nfe_pmd.o nfe_bbmd.o nfe_stsm.o

XRAY_OBJECTS = \
   xray_fourier.o xray_globals.o xray_interface.o xray_utils.o \
   xray_reciprocal_space.o assert.o wallclock.o bulk_solvent.o

ifneq ("$(wildcard ~/.amber_extra_flags.inc)", "")
  include ~/.amber_extra_flags.inc
endif

install: configured_serial $(BINDIR)/pmemd$(SFX)

parallel: configured_parallel $(BINDIR)/pmemd.MPI$(SFX)
	@echo " "

.PHONY: cuda_serial cuda_SPFP cuda_DPFP cuda_lib_clean
cuda_serial: configured_cuda cuda_SPFP cuda_DPFP

cuda_SPFP:
	$(MAKE) cuda_lib_clean
	$(MAKE) $(BINDIR)/pmemd.cuda_SPFP
	cd $(BINDIR) && ln -f -s pmemd.cuda_SPFP$(SFX) pmemd.cuda$(SFX)

cuda_DPFP:
	$(MAKE) cuda_lib_clean
	$(MAKE) $(BINDIR)/pmemd.cuda_DPFP$(SFX)

.PHONY: cuda_parallel
cuda_parallel: configured_cuda configured_parallel cuda_parallel_SPFP cuda_parallel_DPFP

cuda_parallel_SPFP:
	$(MAKE) cuda_lib_clean
	$(MAKE) $(BINDIR)/pmemd.cuda_SPFP.MPI$(SFX)
	cd $(BINDIR) && ln -f -s pmemd.cuda_SPFP.MPI$(SFX) pmemd.cuda.MPI$(SFX);\

cuda_parallel_DPFP:
	$(MAKE) cuda_lib_clean
	$(MAKE) $(BINDIR)/pmemd.cuda_DPFP.MPI$(SFX)

$(BINDIR)/pmemd$(SFX): $(OBJS) $(EMIL) $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) \
    $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) $(XRAY_OBJECTS) -L$(LIBDIR) \
    $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) $(NFE_OBJECTS) $(PLUMED_LOAD)

$(BINDIR)/pmemd.MPI$(SFX): $(OBJS) $(EMIL) $(NFE_OBJECTS) \
    $(PLUMED_DEPENDENCIES) $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) $(XRAY_OBJECTS) -L$(LIBDIR) \
    $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) $(NFE_OBJECTS) $(PLUMED_LOAD)

$(BINDIR)/pmemd.cuda_SPFP$(SFX): $(OBJS) cuda_spfp_libs $(EMIL) \
    $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
    $(XRAY_OBJECTS) \
    $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) \
    $(NFE_OBJECTS) $(PLUMED_LOAD)

$(BINDIR)/pmemd.cuda_SPFP.MPI$(SFX): $(OBJS) cuda_spfp_libs $(EMIL) \
    $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
    $(XRAY_OBJECTS) \
    $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) \
    $(NFE_OBJECTS) $(PLUMED_LOAD)

$(BINDIR)/pmemd.cuda_DPFP$(SFX): $(OBJS) cuda_dpfp_libs $(EMIL) \
    $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
    $(XRAY_OBJECTS) \
    $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) \
    $(NFE_OBJECTS) $(PLUMED_LOAD)

$(BINDIR)/pmemd.cuda_DPFP.MPI$(SFX): $(OBJS) cuda_dpfp_libs $(EMIL) \
    $(NFE_OBJECTS) $(PLUMED_DEPENDENCIES) $(XRAY_OBJECTS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
    $(XRAY_OBJECTS) \
    $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIBF) $(LDFLAGS) $(PMEMD_FLIBSF) \
    $(NFE_OBJECTS) $(PLUMED_LOAD)


# PLUMED
ifndef PLUMED_INCLUDE_FILE
Plumed.o: Plumed.c
	$(PMEMD_CC) $(CFLAGS) $(AMBERCFLAGS) -o Plumed.o -c Plumed.c
endif

boost-gamma-distribution.o: boost-gamma-distribution.cpp
	$(CXX) $(CXXFLAGS) $(AMBERCFLAGS) -I"$(BASEDIR)/AmberTools/src/boost" -o $@ -c $^

EMIL:
	$(MAKE) -C ../../../AmberTools/src/emil install

depends:
	../makef90depends > f90depends

clean:
	$(RM) -f *.o *.obj *.mod pmemd.MPI$(SFX) pmemd$(SFX) pmemd.cuda_SPFP$(SFX) pmemd.cuda_SPFP.MPI$(SFX) pmemd.cuda_DPFP$(SFX) pmemd.cuda_DPFP.MPI$(SFX) pmemd.cuda$(SFX) pmemd.cuda.MPI$(SFX) *.d work.pc*
	$(MAKE) -C ./cuda clean

cuda_lib_clean::
	$(MAKE) -C ./cuda clean

# Control the suffixes available; this essentially eliminates built-in
# inference rules, done here to improve portability.
.SUFFIXES:
.SUFFIXES: .F90 .c $(OSFX)

.F90$(OSFX):
	@echo "[PMEMD]  FC $<"
	$(VB)$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.F90

.c$(OSFX):
	@echo "[PMEMD]  CC $<"
	$(VB)$(PMEMD_CC) $(PMEMD_COPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.c

# compile pmemd prmtop_dat at lower optimization for buggy gnu 5.x: see bug 303.
prmtop_dat$(OSFX): prmtop_dat.F90
	@echo "[PMEMD]  FC $<"
	$(VB)$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) \
	    $(PMEMD_GNU_BUG303) -c prmtop_dat.F90

# compile pme_setup at lower optimization for apparent buggy intel compilers -
# at least 18.0.x: see bug 360.
pme_setup$(OSFX): pme_setup.F90
	@echo "[PMEMD]  FC $<"
	$(VB)$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_INTEL_BUG360) \
	    $(PMEMD_CU_DEFINES) $(NETCDFINC) -c pme_setup.F90

cuda_spfp_libs:
	$(MAKE) -C ./cuda -e PREC_MODEL=-Duse_SPFP cuda.a

cuda_dpfp_libs:
	$(MAKE) -C ./cuda -e PREC_MODEL=-Duse_DPFP cuda.a

AmberNetcdf$(OSFX): ../../../AmberTools/src/lib/AmberNetcdf.F90
	$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c -o AmberNetcdf$(OSFX) ../../../AmberTools/src/lib/AmberNetcdf.F90

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE)!" ;\
	    exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type parallel!" ;\
		echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
		exit 2 ;\
	} ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -mpi." ;\
		exit 2 ;\
	fi ;\
	)
	@(if grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type cuda, not serial!" ;\
		echo "  Rerun $(CONFIG_COMMAND) without -cuda." ;\
		exit 2 ;\
	fi ;\
	)

configured_cuda: configured
	@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type cuda!" ;\
		echo "  Rerun $(CONFIG_COMMAND) and specify a cuda implementation." ;\
		exit 2 ;\
		} ;\
	)

## Dependencies for f90 modules and include files generated by makef90depends:

include f90depends
