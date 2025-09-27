SHELL=/bin/sh
AMBER_VERSION=24

include ../config.h

#Cuda Options
#------------
#Precision model to test (options are SPFP, SPDP, DPDP) - default is DPFP
PREC_MODEL=DPFP
#------------

######    Top-level targets, for general use:    ############################

test: clean
	./test_amber_serial.sh

test.parallel: clean
	./test_amber_parallel.sh

test.cuda.serial: clean
	./test_amber_cuda_serial.sh

test.cuda.serial.SPFP: clean
	./test_amber_cuda_serial.sh SPFP

test.cuda.parallel: clean
	./test_amber_cuda_parallel.sh

test.cuda.parallel.SPFP: clean
	./test_amber_cuda_parallel.sh SPFP

test.hip.serial: clean
	./test_amber_hip_serial.sh

test.hip.serial.SPFP: clean
	./test_amber_hip_serial.sh SPFP

test.hip.parallel: clean
	./test_amber_hip_parallel.sh

test.hip.parallel.SPFP: clean
	./test_amber_hip_parallel.sh SPFP

######    Intermediate-level targets for Amber:    ##################

test.serial.pmemd: test.serial.pmemd.basic test.serial.pmemd.TI \
	test.serial.emil.pmemd test.serial.pmemd.nfe  \
	test.serial.pmemd.amd test.serial.pmemd.kmmd  \
	test.serial.pmemd.middle-scheme

test.parallel.pmemd: test.parallel.pmemd.basic test.parallel.pmemd.TI \
	test.pmemd.REM test.parallel.pmemd.nfe test.pmemd.neb \
	test.parallel.pmemd.amd test.parallel.pmemd.mbar test.parallel.pmemd.kmmd \
	test.parallel.pmemd.middle-scheme

test.cuda.serial2: test.pmemd.cuda

test.parallel.cuda: test.pmemd.cuda.MPI 

test.hip.serial2: test.pmemd.hip

test.parallel.hip: test.pmemd.hip.MPI

######   More specialized rules; not for general use ################

test.parallel.4procs:
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_gb_2rep && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_wat && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_wat_pv && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd hrem_wat_pv && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/TempRem && ./Run.cnste_remd

test.parallel.8procs:
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_gb_4rep && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd h_rem && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd multid_remd && ./Run.multirem
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/trpcage_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/gact_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/gact_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/pHREMrst && ./Run.pHremdrst
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/EREM && ./Run.Eremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/EREMrst && ./Run.Eremdrst
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/Explicit_EREM && ./Run.Eremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phtremd/explicit && ./Run.phtremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phtremd/explicitRst && ./Run.phtremdRst

test.parallel.16procs:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd ephtremd/explicit && ./Run.ephtremd
	-export TESTsander=$(BINDIR)/pmemd.MPI; cd hrem_random_partner && ./Run.hremd

test.pmemd.REM:
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_gb_2rep && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_wat && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_wat_pv && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd hrem_wat_pv && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd rem_gb_4rep && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd h_rem && ./Run.rem
	-export TESTsander='$(BINDIR)/pmemd.MPI' && cd multid_remd && ./Run.multirem
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/trpcage_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/gact_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rxsgld/gact_4rep && ./Run.rxsgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/pHREM && ./Run.pHremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/pHREMrst && ./Run.pHremdrst
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/TempRem && ./Run.cnstph_remd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph_remd/Explicit_pHREM && ./Run.pHremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/EREM && ./Run.Eremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/EREMrst && ./Run.Eremdrst
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/TempRem && ./Run.cnste_remd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste_remd/Explicit_EREM && ./Run.Eremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phtremd/implicit && ./Run.phtremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phtremd/explicit && ./Run.phtremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phtremd/explicitRst && ./Run.phtremdRst
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd ephtremd/implicit && ./Run.ephtremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd ephtremd/explicit && ./Run.ephtremd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phmd_remd/implicit/ && ./Run.phmdremd
	-export TESTsander=$(BINDIR)/pmemd.MPI; cd hrem_random_partner && ./Run.hremd

test.serial.pmemd.middle-scheme:
	-export TESTsander='$(BINDIR)/pmemd'; cd middle-scheme/4096wat && ./Run.pure_wat_middle_scheme
	-export TESTsander='$(BINDIR)/pmemd'; cd middle-scheme/DNA7 && ./Run.DNA_middle_scheme
	-export TESTsander='$(BINDIR)/pmemd'; cd middle-scheme/ETH && ./Run.ETH_middle_scheme

test.serial.pmemd.amd:
	-export TESTsander='$(BINDIR)/pmemd'; cd gamd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='$(BINDIR)/pmemd'; cd amd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='$(BINDIR)/pmemd'; cd scaledMD && $(MAKE) -k test OPT=pmemd

test.serial.pmemd.basic:
	-export TESTsander='$(BINDIR)/pmemd'; cd nonper && ./Run.egb6
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.pure_wat
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.vrand
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.frcdmp
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat_oct && ./Run.pure_wat_oct
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.pure_wat_fswitch
	-export TESTsander='$(BINDIR)/pmemd'; cd 4096wat && ./Run.pure_wat_efield
	-export TESTsander='$(BINDIR)/pmemd'; cd alp && ./Run.alp
	-export TESTsander='$(BINDIR)/pmemd'; cd cytosine && ./Run.cytosine
	-export TESTsander='$(BINDIR)/pmemd'; cd dhfr && ./Run.dhfr
	-export TESTsander='$(BINDIR)/pmemd'; cd dhfr && ./Run.dhfr.min
	-export TESTsander='$(BINDIR)/pmemd'; cd dhfr && ./Run.dhfr.noshake
	-export TESTsander='$(BINDIR)/pmemd'; cd dhfr && ./Run.dhfr.bussi
	-export TESTsander='$(BINDIR)/pmemd'; cd dhfr && ./Run.dhfr_reweight
	-export TESTsander='$(BINDIR)/pmemd'; cd ff14ipq && ./Run.ff14ipq
	-export TESTsander='$(BINDIR)/pmemd'; cd variable_14 && ./Run.variable_14_ntb1
	-export TESTsander='$(BINDIR)/pmemd'; cd variable_14 && ./Run.variable_14_igb1
	-export TESTsander='$(BINDIR)/pmemd'; cd trx && ./Run.trx
	-export TESTsander='$(BINDIR)/pmemd'; cd trx && ./Run.trx.cpln.pmemd
	-export TESTsander='$(BINDIR)/pmemd'; cd gb_rna && ./Run.gbrna
	-export TESTsander='$(BINDIR)/pmemd'; cd gb_rna && ./Run.gbrna.min
	-export TESTsander='$(BINDIR)/pmemd'; cd gb_rna && ./Run.gbrna.ln
	-export TESTsander='$(BINDIR)/pmemd'; cd circ_dna && ./Run.circdna
	-export TESTsander='$(BINDIR)/pmemd'; cd gb2_trx && ./Run.trxox.nogbsa && ./Run.trxox.nogbsa prmtop_an
	-export TESTsander='$(BINDIR)/pmemd'; cd gb7_trx && ./Run.trxox_md
	-export TESTsander='$(BINDIR)/pmemd'; cd gb8_trx && ./Run.trxox
	-export TESTsander='$(BINDIR)/pmemd'; cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-export TESTsander='$(BINDIR)/pmemd'; cd mcwat && ./Run.mcwat
	-export TESTsander='$(BINDIR)/pmemd'; cd RAMD && ./Run.ramd
	-export TESTsander='$(BINDIR)/pmemd'; cd ff19SB && ./Run.pmemd
# gbneck2nu test. TODO: make better alignment for those lines
	-export TESTsander='$(BINDIR)/pmemd'; cd gbneck2nu/1hji && ./Run.1hji_md
	-export TESTsander='$(BINDIR)/pmemd'; cd gbneck2nu/modified_nu/3p4a/ && ./Run.3p4a
# end gbneck2nu test
	-export TESTsander='$(BINDIR)/pmemd'; cd alpb_trx && ./Run.trxox.nogbsa
#	-export TESTsander='$(BINDIR)/pmemd'; cd gb1_cox2 && ./Run.cox2
	-export TESTsander='$(BINDIR)/pmemd'; cd gbsa_xfin && ./Run.gbsa1
	-export TESTsander='$(BINDIR)/pmemd'; cd gbsa_xfin && ./Run.gbsa3
	-export TESTsander='$(BINDIR)/pmemd'; cd tip4p && ./Run.tip4p
	-export TESTsander='$(BINDIR)/pmemd'; cd tip4p && ./Run.tip4p_mcbar
	-export TESTsander='$(BINDIR)/pmemd'; cd tip4p && ./Run.tip4p_nve
	-export TESTsander='$(BINDIR)/pmemd'; cd tip5p && ./Run.tip5p
	-export TESTsander='$(BINDIR)/pmemd'; cd tip5p && ./Run.tip5p_nve
# constant pH test
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstph/implicit && ./Run.cnstph
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstph/explicit && ./Run.cnstph
	-export TESTsander='$(BINDIR)/pmemd'; cd phmd/implicit && ./Run.phmd
# constant E test
	-export TESTsander='$(BINDIR)/pmemd'; cd cnste/implicit && ./Run.cnste
	-export TESTsander='$(BINDIR)/pmemd'; cd cnste/explicit && ./Run.cnste
# constant pH,E test
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstphe/implicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstphe/explicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstphe_merged/implicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd'; cd cnstphe_merged/explicit && ./Run.cnstphe
#CHARMM TESTS
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-export TESTsander='$(BINDIR)/pmemd'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md
#IPS TESTS
	-export TESTsander='$(BINDIR)/pmemd'; cd gact_ips && ./Run.ips
# Constant surface tension (csurften > 0)
	-export TESTsander='$(BINDIR)/pmemd'; cd csurften && ./Run.csurften_z-dir
# Lennard-Jones 12-6-4 potential
	-export TESTsander='$(BINDIR)/pmemd'; cd lj_12-6-4 && ./Run.12-6-4
# Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-export TESTsander='$(BINDIR)/pmemd'; cd csurften && ./Run.csurften_z-dir_npt_3
#    NMR Opt Test Cases
	-export TESTsander='$(BINDIR)/pmemd'; cd nmropt && $(MAKE) pmemd_compat
	-export TESTsander='$(BINDIR)/pmemd'; cd rdc && ./Run.nmr
# NetCDF Traj/Restart tests
	-export TESTsander='$(BINDIR)/pmemd'; cd netcdf && $(MAKE) -k test OPT=$(NETCDF)
#    Directional pressure scaling (baroscalingdir > 0)
	-export TESTsander='$(BINDIR)/pmemd'; cd baroscalingdir && ./Run.baroscalingdir_z-dir
#    Xray
	@(if [ -n "$(MKL)" ]; then \
          $(MAKE) test.serial.pmemd.xray ;\
	fi ;\
	)
# Velocity write
	-export TESTsander='$(BINDIR)/pmemd'; cd velocity_write && ./Run.vel1
	-export TESTsander='$(BINDIR)/pmemd'; cd velocity_write && ./Run.vel2
# Custom virtual sites
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/tip4p && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/tip4p && ./Run.npt
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/tip5p && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/tip5p && ./Run.npt
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/Methanol && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/Methanol && ./Run.npt
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/BromoBenzene && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/BromoBenzene && ./Run.npt
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/DimethylEther && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/DimethylEther && ./Run.npt
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/BenzenePi && ./Run.ec
	-export TESTsander='$(BINDIR)/pmemd'; cd virtual_sites/BenzenePi && ./Run.npt

test.serial.pmemd.xray:
	-export TESTsander='$(BINDIR)/pmemd'; cd xray/2wcu && ./Run.2wcu
	-export TESTsander='$(BINDIR)/pmemd'; cd xray/480d && ./Run.480d

test.serial.pmemd.TI:
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/gasphase && ./Run.gasphase
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/campTI && ./Run.campTI
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/pheMTI && ./Run.0
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/pheMTI && ./Run.1
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/pheMTI && ./Run.lambda0
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/pheMTI && ./Run.lambda1
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/sodium && ./Run.sodium
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/ti_ggcc && ./Run.test1
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/ti_ggcc && ./Run.test2
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/softcore && ./Run_sc
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/softcore/tishake2vsite && ./Run.tishake2vsite
	-export TESTsander='$(BINDIR)/pmemd'; cd pmemdTI/softcore/tishake2vsite && ./Run.tishake2vsite_vacuum

#2 million atom tests are not run by default because they need >3GB Memory
test.serial.big:
	-cd 2_million_atom && ./Run.2mil.min
	-cd 2_million_atom && ./Run.2mil.md
	-cd 2_million_atom && ./Run.2mil.rst.md

test.serial.pmemd.big:
	-export TESTsander='$(BINDIR)/pmemd'; cd 2_million_atom && ./Run.2mil.min
	-export TESTsander='$(BINDIR)/pmemd'; cd 2_million_atom && ./Run.2mil.md
	-export TESTsander='$(BINDIR)/pmemd'; cd 2_million_atom && ./Run.2mil.rst.md


test.serial.pmemd.sgld:
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgld
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgldfp
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgldg
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgmdg
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgldbd2
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgldbd3
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgldbd4
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgmdbd2
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgmdbd3
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gact && ./Run.sgmdbd4
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgld
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgmdg
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgldbd3
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgldbd4
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgmdbd3
	-export TESTsander='$(BINDIR)/pmemd'; cd sgld/gb_rna && ./Run.gbrna.sgmdbd4

test.parallel.pmemd.sgld:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgldfp
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgldg
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgmdg
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgldbd2
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgldbd3
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgldbd4
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgmdbd2
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgmdbd3
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gact && ./Run.sgmdbd4
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgld
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgldfp
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgldg
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgmdg
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgldbd3
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgldbd4
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgmdbd3
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd sgld/gb_rna && ./Run.gbrna.sgmdbd4

test.parallel.pmemd.big:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 2_million_atom && ./Run.2mil.rst.md

test.parallel.pmemd.middle-scheme:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd middle-scheme/4096wat && ./Run.pure_wat_middle_scheme
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd middle-scheme/4096wat && ./Run.pure_wat_middle_scheme_midpoint
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd middle-scheme/DNA7 && ./Run.DNA_middle_scheme
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd middle-scheme/ETH && ./Run.ETH_middle_scheme
	#-export TESTsander='$(BINDIR)/pmemd.MPI'; cd middle-scheme/ETH && ./Run.ETH_middle_scheme_midpoint

test.parallel.pmemd.amd:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gamd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd amd && $(MAKE) -k test OPT=pmemd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd scaledMD && $(MAKE) -k test OPT=pmemd

test.parallel.pmemd.mbar:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd mbar_MPI && ./Run.sh

test.parallel.pmemd.basic:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd nonper && ./Run.egb6
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.pure_wat
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.pure_wat_fswitch
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.pure_wat_efield
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.pure_wat_nmr_temp_reg
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.vrand
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.frcdmp
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat_oct && ./Run.pure_wat_oct
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd alp && ./Run.alp
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cytosine && ./Run.cytosine
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd dhfr && ./Run.dhfr
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd dhfr && ./Run.dhfr.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd dhfr && ./Run.dhfr.noshake
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd ff14ipq && ./Run.ff14ipq
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd variable_14 && ./Run.variable_14_ntb1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd trx && ./Run.trx
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd trx && ./Run.trx.cpln.pmemd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb_rna && ./Run.gbrna
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb_rna && ./Run.gbrna.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb_rna && ./Run.gbrna.ln
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd circ_dna && ./Run.circdna
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb2_trx && ./Run.trxox.nogbsa
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb7_trx && ./Run.trxox_md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb8_trx && ./Run.trxox
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb8_trx && ./Run.trxox_md && ./Run.trxox_md prmtop_an
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gbneck2nu/1hji && ./Run.1hji_md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd alpb_trx && ./Run.trxox.nogbsa
#	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gb1_cox2 && ./Run.cox2
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gbsa_xfin && ./Run.gbsa1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gbsa_xfin && ./Run.gbsa3
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd tip4p && ./Run.tip4p
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd tip4p && ./Run.tip4p_mcbar
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd tip4p && ./Run.tip4p_nve
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd tip5p && ./Run.tip5p
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd tip5p && ./Run.tip5p_nve
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd ff19SB && ./Run.pmemd
# constant pH test
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph/implicit && ./Run.cnstph
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstph/explicit && ./Run.cnstph
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd phmd/implicit && ./Run.phmd
# constant E test
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste/implicit && ./Run.cnste
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnste/explicit && ./Run.cnste
# constant pH,E test
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstphe/implicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstphe/explicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstphe_merged/implicit && ./Run.cnstphe
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd cnstphe_merged/explicit && ./Run.cnstphe
#CHARMM TESTS
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr && ./Run.dhfr_charmm.md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap && ./Run.dhfr_charmm.md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm_noshake.md
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.min
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd chamber/md_engine/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_charmm.md
#IPS TESTS
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd gact_ips && ./Run.ips
# Constant surface tension (csurften > 0)
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd csurften && ./Run.csurften_z-dir
# Lennard-Jones 12-6-4 potential
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd lj_12-6-4 && ./Run.12-6-4
# Semi-isotropic pressure coupling (csurften > 0 AND ntp = 3)
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd csurften && ./Run.csurften_z-dir_npt_3
#    NMR Opt Test Cases
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd nmropt && $(MAKE) pmemd_compat
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd rdc && ./Run.nmr
# NetCDF Traj/Restart tests
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd netcdf && $(MAKE) -k test OPT=$(NETCDF)
# MIDPOINT TEST CASES
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd 4096wat && ./Run.pure_wat_midpoint
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd dhfr && ./Run.dhfr_midpoint
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd dhfr && ./Run.dhfr.noshake_midpoint
# Velocity write
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd velocity_write && ./Run.vel1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd velocity_write && ./Run.vel2

#EMIL TI TESTS
test.parallel.emil:
	-export TESTsander='$(BINDIR)/sander.MPI'; cd emil/emil_sander_gbsa && ./Run.emil $(EMIL)

test.parallel.emil.pmemd:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd emil/emil_pmemd_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd emil/emil_pmemd_tip3p && ./Run.emil $(EMIL)

test.serial.emil.pmemd:
	-export TESTsander='$(BINDIR)/pmemd';  unset DO_PARALLEL; cd emil/emil_pmemd_gbsa && ./Run.emil $(EMIL)
	-export TESTsander='$(BINDIR)/pmemd';  unset DO_PARALLEL; cd emil/emil_pmemd_tip3p && ./Run.emil $(EMIL)

#KMMD TESTS
test.serial.pmemd.kmmd:
	-export TESTsander='$(BINDIR)/pmemd'; unset DO_PARALLEL; cd kmmd/kmmd_pmemd    && ./Run.kmmd
	-export TESTsander='$(BINDIR)/pmemd'; unset DO_PARALLEL; cd kmmd/kmmd_pmemd_gb && ./Run.kmmd

test.parallel.pmemd.kmmd:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd kmmd/kmmd_pmemd    && ./Run.kmmd
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd kmmd/kmmd_pmemd_gb && ./Run.kmmd

#TI TESTS
test.parallel.pmemd.TI:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/campTI && ./Run.campTI
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.0
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.lambda0
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/pheMTI && ./Run.lambda1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/sodium && ./Run.sodium
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/ti_ggcc && ./Run.test1
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/ti_ggcc && ./Run.test2
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd pmemdTI/softcore && ./Run_sc

test.serial.pmemd.emap:
	-export TESTsander='$(BINDIR)/pmemd'; cd emap/ && ./Run.emap
	-export TESTsander='$(BINDIR)/pmemd'; cd emap/ && ./Run.emap_boundary

test.parallel.pmemd.emap:
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd emap/ && ./Run.emap
	-export TESTsander='$(BINDIR)/pmemd.MPI'; cd emap/ && ./Run.emap_boundary

test.pmemd.midpoint:
	cd midpoint && $(MAKE) -k test.pmemd.midpoint

test.pmemd.cuda:
	cd cuda && $(MAKE) -k test.pmemd.cuda PREC_MODEL=$(PREC_MODEL)

test.pmemd.cuda.MPI:
	cd cuda && $(MAKE) -k test.pmemd.cuda.MPI PREC_MODEL=$(PREC_MODEL)

test.pmemd.hip:
	cd hip && $(MAKE) -k test.pmemd.hip PREC_MODEL=$(PREC_MODEL)

test.pmemd.hip.MPI:
	cd hip && $(MAKE) -k test.pmemd.hip.MPI PREC_MODEL=$(PREC_MODEL)

test.pmemd.cuda.sgld:
	cd cuda && $(MAKE) -k test.pmemd.cuda.sgld PREC_MODEL=$(PREC_MODEL)

test.pmemd.cuda.MPI.sgld:
	cd cuda && $(MAKE) -k test.pmemd.cuda.sgld PREC_MODEL=$(PREC_MODEL)

test.pmemd.cuda.MPI.rxsgld:
	cd cuda && $(MAKE) -k test.pmemd.cuda.rxsgld PREC_MODEL=$(PREC_MODEL)

test.serial.pmemd.nfe:
	-cd nfe && ./run-pmemd-serial.sh

test.parallel.pmemd.nfe:
	-cd nfe && ./run-pmemd-parallel.sh

test.netcdf:
	-(cd ../AmberTools/src/netcdf-4.3.0 && $(MAKE) check)
	-cd trajene_netcdf && ./Run.trajene $(NETCDF)
	-cd netcdf && $(MAKE) -k test OPT=$(NETCDF)

test.pmemd.neb:
	-cd neb-testcases/neb_gb_partial && ./Run.neb_gb_partial_pmemd
	-cd neb-testcases/neb_gb_full    && ./Run.neb_gb_full_pmemd
	-cd neb-testcases/neb_explicit   && ./Run.neb_explicit_pmemd

#===========================================================================

finished.serial:
	@echo ""
	@echo "Finished serial test suite for Amber ${AMBER_VERSION} at `date`."
	@echo ""

finished.parallel:
	@echo ""
	@echo "Finished parallel test suite for Amber ${AMBER_VERSION} at `date`."
	@echo "Tests ran with DO_PARALLEL=\"${DO_PARALLEL}\"."
	@echo "Some tests require 4 or 8 threads to run, while some will not"
	@echo "run with more than 2.  Please run further parallel tests with the"
	@echo "appropriate number of processors. See $(AMBERHOME)/test/README."
	@echo ""

finished.cuda.serial:
	@echo ""
	@echo "Finished CUDA test suite for Amber ${AMBER_VERSION} at `date`."
	@echo ""

finished.hip.serial:
	@echo ""
	@echo "Finished HIP test suite for Amber ${AMBER_VERSION} at `date`."
	@echo ""

clean:
	-./test_amber_clean.sh

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

#===========================================================================
