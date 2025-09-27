# CUDA REMD test Makefile

test:
	./test_cuda_remd.sh

test.complete: rem_2rep_gb rem_4rep_gb rem_2rep_pme hrem_4rep_gb mrem_4rep_gb test.wat

test.2rep: rem_2rep_gb rem_2rep_pme rem_wat_2 rem_wat_pv hrem_wat_pv test.wat rrem_2rep_gb_boltzmann rrem_2rep_gb_nonboltzmann rem_2rep_hybridsolvent

test.4rep: rem_4rep_gb hrem_4rep_gb mrem_4rep_gb rrem_4rep_gb_boltzmann rrem_4rep_gb_nonboltzmann rem_4rep_hybridsolvent

test.wat: rem_wat_2 rem_wat_pv hrem_wat_pv

rem_2rep_gb::
	@-cd rem_2rep_gb && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_4rep_gb::
	@-cd rem_4rep_gb && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_2rep_pme::
	@-cd rem_2rep_pme && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

hrem_4rep_gb::
	@-cd hrem_4rep_gb && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

mrem_4rep_gb::
	@-cd mrem_4rep_gb && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_wat_2::
	@-cd rem_wat_2 && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_wat_pv::
	@-cd rem_wat_pv && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

hrem_wat_pv::
	@-cd hrem_wat_pv && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rrem_2rep_gb_boltzmann::
	@-cd rrem_2rep_gb_boltzmann && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rrem_4rep_gb_boltzmann::
	@-cd rrem_4rep_gb_boltzmann && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rrem_2rep_gb_nonboltzmann::
	@-cd rrem_2rep_gb_nonboltzmann && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rrem_4rep_gb_nonboltzmann::
	@-cd rrem_4rep_gb_nonboltzmann && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_2rep_hybridsolvent::
	@-cd rem_2rep_hybridsolvent && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

rem_4rep_hybridsolvent::
	@-cd rem_4rep_hybridsolvent && ./Run.rem.sh $(PREC_MODEL) $(NETCDF)

clean:
	@$(MAKE) test.complete PREC_MODEL=clean
