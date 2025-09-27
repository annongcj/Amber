include $(AMBERHOME)/config.h

SHELL=/bin/sh

#Precision model to test (options are SPFP, DPFP) - default is SPFP
#PREC_MODEL=DPFP

######    Top-level targets, for general use:    ############################

test.pmemd.cuda: test.pmemd.cuda.gb.serial test.pmemd.cuda.gb \
		 test.pmemd.cuda.pme test.pmemd.cuda.pme.serial \
		 test.pmemd.cuda.gb.gamd test.pmemd.cuda.pme.gamd test.pmemd.cuda.pme.sgamd\
		 test.pmemd.cuda.VirtualSites test.pmemd.cuda.sgld \
		 test.pmemd.cuda.$(FEP_MODE)

test.pmemd.cuda.BIG: test.pmemd.cuda.pme.big

test.pmemd.cuda.MPI: test.pmemd.cuda.gb test.pmemd.cuda.pme \
		     test.pmemd.cuda.remd test.pmemd.cuda.nfe.parallel test.pmemd.cuda.neb \
		     test.pmemd.cuda.gb.gamd test.pmemd.cuda.pme.gamd \
                   test.pmemd.cuda.gti.remd test.pmemd.cuda.sgld test.pmemd.cuda.rxsgld 

test.pmemd.cuda.gti: test.pmemd.cuda.gti.TI \
		 test.pmemd.cuda.gti.TI_restraint \
		 test.pmemd.cuda.gti.1264 \
		 test.pmemd.cuda.gti.BAT \
		 test.pmemd.cuda.gti.LambdaScheduling \
		 test.pmemd.cuda.gti.SC_Correction \
		 test.pmemd.cuda.gti.NB_EXP \
		 test.pmemd.cuda.gti.SC_Correction \
		 test.pmemd.cuda.gti.others


test.pmemd.cuda.gb.serial:
#Begin Serial Implicit Solvent Tests
	@echo "------------------------------------"
	@echo "Running CUDA Implicit solvent tests."
	@echo "  Precision Model = $(PREC_MODEL)"
	@echo "------------------------------------"
	-cd gb_ala3/ && ./Run.igb1_ntc1_min  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb1_ntc1  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest0_ntt0_igb1_ntc1_hotstart  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest0_ntt0_igb1_ntc1_coldstart  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb1_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb2_ntc1  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb2_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc1  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb7_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb8_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc2_ntr1  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc2_saltcon  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc2_rgbmax  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt0_igb5_ntc2_alpb  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt1_igb1_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt2_igb1_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd gb_ala3/ && ./Run.irest1_ntt3_igb1_ntc2  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/multi-angle/ && ./Run.nmropt_2angle_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/angle_torsion_com/ && ./Run.angdih_gb $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr/ && ./Run.dhfr_charmm.min  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap/ && ./Run.dhfr_charmm.min  $(PREC_MODEL) $(NETCDF)
	-cd nucleosome/ && ./Run_min.1  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.gb:
#Begin Implicit Solvent Tests
	@echo "---------------------------------------------"
	@echo "Running Extended CUDA Implicit solvent tests."
	@echo "      Precision Model = $(PREC_MODEL)"
	@echo "---------------------------------------------"
	-cd trpcage/ && ./Run_md_trpcage  $(PREC_MODEL) $(NETCDF)
	-cd myoglobin/ && ./Run_md_myoglobin  $(PREC_MODEL) $(NETCDF)
	-cd myoglobin/ && ./Run_md_myoglobin_igb7  $(PREC_MODEL) $(NETCDF)
	-cd myoglobin/ && ./Run_md_myoglobin_igb8  $(PREC_MODEL) $(NETCDF)
	-cd myoglobin/ && ./Run_md_myoglobin_igb8_gbsa  $(PREC_MODEL) $(NETCDF)
	-cd myoglobin/ && ./Run_md_myoglobin_igb8_gbsa3  $(PREC_MODEL) $(NETCDF)
	-cd gbsa_xfin/ && ./Run.gbsa3  $(PREC_MODEL) $(NETCDF)
	-cd cnstph/implicit/ && ./Run.cnstph  $(PREC_MODEL) $(NETCDF)
	-cd cnste/implicit/ && ./Run.cnste  $(PREC_MODEL) $(NETCDF)
	-cd cnstphe/implicit/ && ./Run.cnstphe  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr/ && ./Run.dhfr_charmm.md  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap/ && ./Run.dhfr_charmm.md  $(PREC_MODEL) $(NETCDF)
	-cd nucleosome/ && ./Run_md.1  $(PREC_MODEL) $(NETCDF)
	-cd nucleosome/ && ./Run_md.2  $(PREC_MODEL) $(NETCDF)
	-cd amd/rna_gb && ./Run.gb.amd1  $(PREC_MODEL) $(NETCDF)
	-cd amd/rna_gb && ./Run.gb.amd2  $(PREC_MODEL) $(NETCDF)
	-cd amd/rna_gb && ./Run.gb.amd3  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/angle/ && ./Run.nmropt_1angle_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/distance/ && ./Run.dist_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/distance_COM/ && ./Run.distCOM_gb  $(PREC_MODEL) $(NETCDF)
#	-cd nmropt/gb/distance_COM/ && ./Run.distr6av_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/jar_distance_COM/ && ./Run.jar_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/jar_distance/ && ./Run.jar_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/jar_torsion/ && ./Run.jar_torsion  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/nmropt_1_torsion/ && ./Run.nmropt_1_torsion  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/tautp/ && ./Run.nmropt_1tautp_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/temp/ && ./Run.nmropt_1temp_gb  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/gb/angle_torsion_com/ && ./Run.angdih_gb  $(PREC_MODEL) $(NETCDF)
	-cd kmmd/kmmd_pmemd_gb && ./Run.kmmd $(PREC_MODEL) $(NETCDF)


test.pmemd.cuda.gb.gamd:
#Begin Implicit Solvent Tests with GaMD
	-cd gamd/rna_gb && ./Run.gb.gamd1  $(PREC_MODEL) $(NETCDF)
	-cd gamd/rna_gb && ./Run.gb.gamd2  $(PREC_MODEL) $(NETCDF)
	-cd gamd/rna_gb && ./Run.gb.gamd3  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.sgld:
#Begin SGLD Tests 
	-cd sgld/trpcage && ./Run.sgld  $(PREC_MODEL) $(NETCDF)
	-cd sgld/trpcage && ./Run.sgmd  $(PREC_MODEL) $(NETCDF)
	-cd sgld/gact && ./Run.sgld  $(PREC_MODEL) $(NETCDF)
	-cd sgld/gact && ./Run.sgmd  $(PREC_MODEL) $(NETCDF)

#End Implicit Solvent Tests

test.pmemd.cuda.pme.serial:
# ntt=2 and ntt=3 give different
# answers in parallel due to random number differences.
# Do not run in parallel for now to avoid confusion.
	@echo "------------------------------------"
	@echo "Running CUDA Explicit solvent tests."
	@echo "  Precision Model = $(PREC_MODEL)"
	@echo "------------------------------------"
	-cd 4096wat/ && ./Run.vrand  $(PREC_MODEL) $(NETCDF)
	-cd 4096wat_oct/ && ./Run.pure_wat_oct_NVT_NTT2  $(PREC_MODEL) $(NETCDF)
	-cd 4096wat_oct/ && ./Run.pure_wat_oct_NVT_NTT3  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.ntb2_ntt3  $(PREC_MODEL) $(NETCDF)
	-cd ips/ && ./Run.ips_ntt3  $(PREC_MODEL) $(NETCDF)
	-cd npt_ntr_ntt3_oct/ && ./Run.ntp1  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.min  $(PREC_MODEL) $(NETCDF)
	-cd mcwat && ./Run.mcwat $(PREC_MODEL) $(NETCDF)
	-cd RAMD && ./Run.ramd $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_pbc/ && ./Run.dhfr_pbc_charmm_noshake.min  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap_pbc/ && ./Run.dhfr_cmap_pbc_charmm.min  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap_pbc/ && ./Run.dhfr_cmap_pbc_charmm_noshake.min  $(PREC_MODEL) $(NETCDF)
#RCW Explicit solvent constantph appears to suffer from the same issue.
	-cd nfe/ && ./run-pmemd-serial.sh  $(PREC_MODEL) $(NETCDF)
#Newly implemented CMAP test for ff19SB
	-cd ff19SB/ && ./Run  $(PREC_MODEL) $(NETCDF)
#reweight dhfr not working in parallel as yet
	-cd dhfr/ && ./Run.dhfr_reweight  $(PREC_MODEL) $(NETCDF)
#Middle scheme Tests
	-cd ../middle-scheme/4096wat/ && ./Run.pure_wat_middle_scheme.cuda $(PREC_MODEL) $(NETCDF)
	-cd ../middle-scheme/DNA7/ && ./Run.DNA_middle_scheme.cuda $(PREC_MODEL) $(NETCDF)
	-cd ../middle-scheme/ETH/ && ./Run.ETH_middle_scheme.cuda $(PREC_MODEL) $(NETCDF)
	@(if [ -n "$(MKL)" ]; then \
	  $(MAKE) test.pmemd.cuda.pme.xray.serial ;\
	fi ;\
	)

test.pmemd.cuda.pme.xray.serial:
	-cd xray/480d && ./Run.480d  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.pme:
#Begin Explicit Solvent Tests
	@echo "---------------------------------------------"
	@echo "Running Extended CUDA Explicit solvent tests."
	@echo "      Precision Model = $(PREC_MODEL)"
	@echo "---------------------------------------------"
	-cd 4096wat/ && ./Run.pure_wat  $(PREC_MODEL) $(NETCDF)
	-cd 4096wat/ && ./Run.pure_wat_fswitch $(PREC_MODEL) $(NETCDF)
	-cd 4096wat/ && ./Run.pure_wat_efield $(PREC_MODEL) $(NETCDF)
	-cd 4096wat_oct/ && ./Run.pure_wat_oct_NVE  $(PREC_MODEL) $(NETCDF)
	-cd 4096wat_oct/ && ./Run.pure_wat_oct_NVT_NTT1  $(PREC_MODEL) $(NETCDF)
	-cd 4096wat_oct/ && ./Run.pure_wat_oct_NPT_NTT1  $(PREC_MODEL) $(NETCDF)
	-cd large_solute_count/ && ./Run.ntb2_ntt1  $(PREC_MODEL) $(NETCDF)
	-cd jac/ && ./Run.jac  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.ntr1  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.ntb2  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.ntb2_ntt1  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.ntb2_ntt1_ntr1  $(PREC_MODEL) $(NETCDF)
	-cd dhfr/ && ./Run.dhfr.noshake  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_pbc/ && ./Run.dhfr_pbc_charmm_noshake.md  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap_pbc/ && ./Run.dhfr_cmap_pbc_charmm.md  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap_pbc/ && ./Run.dhfr_cmap_pbc_charmm_NPT.md  $(PREC_MODEL) $(NETCDF)
	-cd chamber/dhfr_cmap_pbc/ && ./Run.dhfr_cmap_pbc_charmm_noshake.md  $(PREC_MODEL) $(NETCDF)
	-cd ips/ && ./Run.ips  $(PREC_MODEL) $(NETCDF)
	-cd amd/dhfr_pme && ./Run.pme.amd1  $(PREC_MODEL) $(NETCDF)
	-cd amd/dhfr_pme && ./Run.pme.amd2  $(PREC_MODEL) $(NETCDF)
	-cd amd/dhfr_pme && ./Run.pme.amd3  $(PREC_MODEL) $(NETCDF)
	-cd amd/gact_ips && ./Run.ips.amd1  $(PREC_MODEL) $(NETCDF)
	-cd amd/gact_ips && ./Run.ips.amd2  $(PREC_MODEL) $(NETCDF)
	-cd amd/gact_ips && ./Run.ips.amd3  $(PREC_MODEL) $(NETCDF)
	-cd scaledMD && ./Run.scaledMD  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/angle/ && ./Run.nmropt_1angle_pbc  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/distance/ && ./Run.dist_pbc  $(PREC_MODEL) $(NETCDF)
#Newly implemented COM support for PME including jar
	-cd nmropt/pme/distance_COM/ && ./Run.distCOM_pbc $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/distance_COM/ && ./Run.multi-distCOM-multi-dist_pbc $(PREC_MODEL) $(NETCDF)
#	-cd nmropt/pme/distance_COM/ && ./Run.multi-distr6av_pbc $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/jar_distance_COM/ && ./Run.jar_pbc  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/jar_distance/ && ./Run.jar_pbc  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/jar_torsion/ && ./Run.jar_torsion  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/nmropt_1_torsion/ && ./Run.nmropt_1_torsion  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/angle_torsion_com/ && ./Run.angdihmultCOM  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/tautp/ && ./Run.nmropt_1tautp_pbc  $(PREC_MODEL) $(NETCDF)
	-cd nmropt/pme/temp/ && ./Run.nmropt_1temp_pbc  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_box_nvt  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_box_npt  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_oct_nvt  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_oct_npt  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_box_npt_mcbar  $(PREC_MODEL) $(NETCDF)
	-cd tip4pew/ && ./Run.tip4pew_oct_npt_mcbar  $(PREC_MODEL) $(NETCDF)
	-cd tip5p/ && ./Run.tip5p_box_nvt  $(PREC_MODEL) $(NETCDF)
	-cd tip5p/ && ./Run.tip5p_box_npt  $(PREC_MODEL) $(NETCDF)
	-cd tip5p/ && ./Run.tip5p_oct_nvt  $(PREC_MODEL) $(NETCDF)
	-cd tip5p/ && ./Run.tip5p_oct_npt  $(PREC_MODEL) $(NETCDF)
	-cd HMR/ && ./Run.HMR	$(PREC_MODEL) $(NETCDF)
	-cd lipid_npt_tests/ && ./Run_npt_isotropic.lipid14	$(PREC_MODEL) $(NETCDF)
	-cd lipid_npt_tests/ && ./Run_npt_anisotropic.lipid14	$(PREC_MODEL) $(NETCDF)
	-cd lipid_npt_tests/ && ./Run_npt_semiisotropic.lipid14	$(PREC_MODEL) $(NETCDF)
	-cd cnstph/explicit/ && ./Run.cnstph  $(PREC_MODEL) $(NETCDF)
	-cd cnste/explicit/ && ./Run.cnste  $(PREC_MODEL) $(NETCDF)
	-cd cnstphe/explicit/ && ./Run.cnstphe  $(PREC_MODEL) $(NETCDF)
	-cd phmd/explicit/ && ./Run.phmd $(PREC_MODEL) $(NETCDF)
	-cd kmmd/kmmd_pmemd/ && ./Run.kmmd $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.pme.gamd:
	-cd gamd/dhfr_pme && ./Run.pme.gamd1  $(PREC_MODEL) $(NETCDF)
	-cd gamd/dhfr_pme && ./Run.pme.gamd2  $(PREC_MODEL) $(NETCDF)
	-cd gamd/dhfr_pme && ./Run.pme.gamd3  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.pme.sgamd:
	-cd gamd/LiGaMD && ./Run.liggamd  $(PREC_MODEL) $(NETCDF)
	-cd gamd/PepGaMD && ./Run.pepgamd $(PREC_MODEL) $(NETCDF)
	-cd gamd/PPIGaMD && ./Run.ppigamd $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.remd:
	-cd remd && $(MAKE) test
	-cd cnstph_remd/implicit && ./Run.pHremd  $(PREC_MODEL) $(NETCDF)
	-cd cnstph_remd/explicit && ./Run.pHremd  $(PREC_MODEL) $(NETCDF)
#	-cd remd/rem_gb_2rep && ./Run.rem  $(PREC_MODEL) $(NETCDF)
#	-cd remd/rem_wat && ./Run.rem  $(PREC_MODEL) $(NETCDF)
#	-cd remd/rem_gb_4rep && ./Run.rem  $(PREC_MODEL) $(NETCDF)
#	-cd remd/multid_remd && ./Run.multirem  $(PREC_MODEL) $(NETCDF)
#	-cd remd/h_rem && ./Run.rem  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.rxsgld:
	-cd rxsgld && $(MAKE) test

test.pmemd.cuda.pme.big:
#Large test cases which have problems on <=3GB in parallel.
	-cd FactorIX_NVE/ && ./Run.FactorIX_NVE  $(PREC_MODEL) $(NETCDF)
	-cd cellulose/ && ./Run.cellulose_nvt_256_128_128  $(PREC_MODEL) $(NETCDF)
	-cd cellulose/ && ./Run.cellulose_nvt  $(PREC_MODEL) $(NETCDF)
	-cd cellulose/ && ./Run.cellulose_npt  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.nfe.parallel:
	-cd nfe/ && ./run-pmemd-parallel.sh  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.neb:
	-cd neb-testcases/neb_gb_partial && ./Run.neb_gb_partial $(PREC_MODEL) $(NETCDF)
	-cd neb-testcases/neb_gb_full    && ./Run.neb_gb_full    $(PREC_MODEL) $(NETCDF)
# DAC, 4/16/2024: skip this test, since the box size is too small for
# pmemd.CUDA
# 	-cd neb-testcases/neb_explicit   && ./Run.neb_explicit   $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.afe:
#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA AFE free energy tests."
	-cd AFE/TI/campTI && ./Run.campTI  $(PREC_MODEL) $(NETCDF)
	-cd AFE/TI/softcore/complex && ./Run.complex  $(PREC_MODEL) $(NETCDF)
	-cd AFE/TI/softcore/electro && ./Run.electro  $(PREC_MODEL) $(NETCDF)
	-cd AFE/TI/softcore/dynlmb && ./Run.dynlmb  $(PREC_MODEL) $(NETCDF)
	-cd AFE/MBAR/campMBAR && ./Run.campMBAR  $(PREC_MODEL) $(NETCDF)
	-cd AFE/MBAR/softcore/complex && ./Run.complex  $(PREC_MODEL) $(NETCDF)
	-cd AFE/MBAR/softcore/electro && ./Run.electro  $(PREC_MODEL) $(NETCDF)

test.pmemd.cuda.gti.TI:
#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI free energy tests."
	-cd gti/Na/ && ./Run.NVE  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.NVT_ntt1  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.NPT_MC  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/Na/ && ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.NVE  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.NVT_ntt1  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.NPT_MC  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/ligand/ && ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.NVE  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.NVT_ntt1  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.NPT_MC  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/complex/ && ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/AFE_TI/campTI && ./Run.campTI  $(PREC_MODEL)
	-cd gti/AFE_TI/softcore/complex && ./Run.complex  $(PREC_MODEL)
	-cd gti/AFE_TI/softcore/electro && ./Run.electro  $(PREC_MODEL)
	-cd gti/AFE_TI/softcore/dynlmb && ./Run.dynlmb  $(PREC_MODEL)
	-cd gti/AFE_MBAR/campMBAR && ./Run.campMBAR  $(PREC_MODEL)  SCH=222 ./Run.SC_NVT_MBAR SPFP
	-cd gti/AFE_MBAR/softcore/complex && ./Run.complex  $(PREC_MODEL)
	-cd gti/AFE_MBAR/softcore/electro && ./Run.electro  $(PREC_MODEL)
	-cd gti/tishake2vsite && Run.tishake2vsite  $(PREC_MODEL)

test.pmemd.cuda.gti.TI_restraint:
#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI NMR Restratin tests."
	-cd gti/NMR_Restraint/SimpleTorsion && CASE=A ./Run  $(PREC_MODEL)
	-cd gti/NMR_Restraint/SimpleTorsion && CASE=B ./Run  $(PREC_MODEL)
	-cd gti/NMR_Restraint/SimpleTorsion && CASE=C ./Run  $(PREC_MODEL)

test.pmemd.cuda.gti.1264:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI 12-6-4 tests."
	-cd gti/1264/MgMD && ./Run.NVT  $(PREC_MODEL)
	-cd gti/1264/2ionMD && ./Run.NVT  $(PREC_MODEL)
	-cd gti/1264/Mg2Zn && ./Run.NVT  $(PREC_MODEL)
	-cd gti/1264/SimpleTI && ./Run  $(PREC_MODEL)

test.pmemd.cuda.gti.BAT:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI BAT tests."
	-cd gti/BAT/Simple && DIR=for ./Run  $(PREC_MODEL)
	-cd gti/BAT/Simple && DIR=back ./Run  $(PREC_MODEL)

test.pmemd.cuda.gti.SC_Correction:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI free energy SC correction tests."
	-cd gti/SC_Correction/methane_2_methanol/ && SC_CORR=0 ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/SC_Correction/methane_2_methanol/ && SC_CORR=-1 ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/SC_Correction/methane_2_methanol/ && SC_CORR=1 ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/SC_Correction/methane_2_methanol/ && SC_CORR=2 ./Run.NVT_ntt3  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=0 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=-1 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=1 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=2 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=0 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=-1 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=1 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=2 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=0 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=-1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/ligand/ && SC_CORR=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=0 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=-1 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=1 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=2 ./Run.SC_NVE  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=0 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=-1 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=1 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=2 ./Run.SC_NVT  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=0 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=-1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/SC_Correction/complex/ && SC_CORR=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)

test.pmemd.cuda.gti.NB_EXP:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI free energy NB exponent tests."
	-cd gti/NB_EXP/ligand/ && NB_N=6 NB_M=6 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=1 NB_M=6 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=2 NB_M=6 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=4 NB_M=6 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=6 NB_M=4 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=1 NB_M=4 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=2 NB_M=4 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=4 NB_M=4 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=6 NB_M=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=1 NB_M=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=2 NB_M=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=4 NB_M=2 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=6 NB_M=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=1 NB_M=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=2 NB_M=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)
	-cd gti/NB_EXP/ligand/ && NB_N=4 NB_M=1 ./Run.SC_NVT_MBAR  $(PREC_MODEL)

test.pmemd.cuda.gti.LambdaScheduling:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI Lambda Scheduling."
	-cd gti/LambdaScheduling/complex/ && SCH=222 ./Run.SC_NVT_MBAR $(PREC_MODEL)
	-cd gti/LambdaScheduling/complex/ && SCH=333 ./Run.SC_NVT_MBAR $(PREC_MODEL)
	-cd gti/LambdaScheduling/complex/ && SCH=444 ./Run.SC_NVT_MBAR $(PREC_MODEL)
	-cd gti/LambdaScheduling/complex/ && SCH=22T ./Run.SC_NVT_MBAR $(PREC_MODEL)

test.pmemd.cuda.gti.remd:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI Lambda Replica exchange Scheduling."
	-cd gti/lambda_remd/multi-window/ && NREP=4 ./Run $(PREC_MODEL)
	-cd gti/lambda_remd/multi-window/ && NREP=5 ./Run $(PREC_MODEL)
	-cd gti/lambda_remd/multi-window/ && NREP=8 ./Run $(PREC_MODEL)
	-cd gti/lambda_remd/multi-window/ && NREP=11 ./Run $(PREC_MODEL)
#	-cd gti/lambda_remd/diphenyltoluene/ && ./Run $(PREC_MODEL)

test.pmemd.cuda.gti.others:
	#Begin GTI Tests
	@echo "------------------------------------"
	@echo "Running CUDA GTI free energy other tests."
	-cd gti/REAF/ligand && LAMBDA=0 TAU=0 ./Run.REAF $(PREC_MODEL)
	-cd gti/REAF/ligand && LAMBDA=0 TAU=0.292893218813 ./Run.REAF $(PREC_MODEL)
	-cd gti/REAF/ligand && LAMBDA=0.5 TAU=0 ./Run.REAF $(PREC_MODEL)
	-cd gti/REAF/ligand && LAMBDA=1 TAU=0.292893218813 ./Run.REAF $(PREC_MODEL)
	-cd gti/REAF/ligand && LAMBDA=1 TAU=0 ./Run.REAF $(PREC_MODEL)
	-cd gti/SHAKE/cyclopentane_2_thf && LAMBDA=0 ./Run $(PREC_MODEL)
	-cd gti/SHAKE/cyclopentane_2_thf && LAMBDA=0.5 ./Run $(PREC_MODEL)
	-cd gti/SHAKE/cyclopentane_2_thf && LAMBDA=1 ./Run $(PREC_MODEL)
	-cd gti/SHAKE/JAK2 && LAMBDA=0 ./Run $(PREC_MODEL)
	-cd gti/SHAKE/JAK2 && LAMBDA=0.5 ./Run $(PREC_MODEL)
	-cd gti/SHAKE/JAK2 && LAMBDA=1.0 ./Run $(PREC_MODEL)
	-cd gti/TargetVolume/ligand && ./Run $(PREC_MODEL)


test.pmemd.cuda.VirtualSites:
	#Begin tests
	@echo "------------------------------------"
	@echo "Running CUDA Virtual Site tests."
	-cd virtual_sites/tip4p && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/tip4p && ./Run.npt $(PREC_MODEL)
	-cd virtual_sites/tip5p && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/tip5p && ./Run.npt $(PREC_MODEL)
	-cd virtual_sites/Methanol && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/Methanol && ./Run.npt $(PREC_MODEL)
	-cd virtual_sites/BromoBenzene && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/BromoBenzene && ./Run.npt $(PREC_MODEL)
	-cd virtual_sites/DimethylEther && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/DimethylEther && ./Run.npt $(PREC_MODEL)
	-cd virtual_sites/BenzenePi && ./Run.ec $(PREC_MODEL)
	-cd virtual_sites/BenzenePi && ./Run.npt $(PREC_MODEL)

#End Explicit Solvent Tests
