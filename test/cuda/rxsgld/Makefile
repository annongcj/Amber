# CUDA REMD test Makefile

test:
	./test_cuda_rxsgld.sh

test.complete: trpcage gact


test.4rep: trpcage  gact


trpcage::
	@-cd trpcage && ./Run.rxsgld.sh $(PREC_MODEL) $(NETCDF)

gact::
	@-cd gact && ./Run.rxsgld.sh $(PREC_MODEL) $(NETCDF)


clean:
	@$(MAKE) test.complete PREC_MODEL=clean
