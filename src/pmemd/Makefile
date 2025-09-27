SHELL=/bin/sh
AMBER=Amber24
CONFIG_FILE=../config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA
MIC_TOKEN=mmic
MIC2_TOKEN=DMIC2
MIC_OFFLOAD_TOKEN=DMIC_offload
HYBRID_TOKEN=DHybrid

include $(CONFIG_FILE)

.NOTPARALLEL: configured_serial configured_parallel configured_cuda configured_mic configured_mic_offload configured_mic2

install: $(INSTALLTYPE)

openmp: serial

cray_serial: serial

openmp:  serial

serial: configured_serial
	$(MAKE) -C src/ install
	@echo "Installation of pmemd complete"

parallel: configured_parallel
	$(MAKE) -C src/ parallel
	@echo "Installation of pmemd.MPI complete"

cuda_serial: configured_cuda
	$(MAKE) -C src/ cuda_serial
	@echo "Installation of pmemd.cuda complete"

cuda_parallel: configured_cuda configured_parallel
	$(MAKE) -C src/ cuda_parallel
	@echo "Installation of pmemd.cuda.MPI complete"

mic: configured_mic
	$(MAKE) -C src/ mic 
	@echo "Installation of pmemd.mic_native complete"

mic_parallel: configured_mic configured_parallel
	$(MAKE) -C src/ mic_parallel
	@echo "Installation of pmemd.mic_native.MPI complete"

mic2: configured_mic2
	$(MAKE) -C src/ mic2
	@echo "Installation of pmemd.MPI complete" 

mic_offload: configured_mic_offload
	$(MAKE) -C src/ mic_offload                    
	@echo "Installation of pmemd.mic_offload.MPI complete" 

hybrid: configured_hybrid
	$(MAKE) -C src/ hybrid
	@echo "Installation of pmemd.hybrid.MPI complete" 

mic_hybrid: configured_hybrid configured_mic
	$(MAKE) -C src/ mic_hybrid
	@echo "Installation of pmemd.hybrid.MPI.mic complete"

uninstall.serial: clean
	-rm -f $(BINDIR)/pmemd$(SFX)

uninstall.cuda: clean
	-rm -f $(BINDIR)/pmemd.cuda$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPFP$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPFP.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPXP$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_SPXP.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_DPFP$(SFX)
	-rm -f $(BINDIR)/pmemd.cuda_DPFP.MPI$(SFX)

uninstall.mic: clean
	-rm -f $(BINDIR)/pmemd.mic_native$(SFX)
	-rm -f $(BINDIR)/pmemd.mic_native.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.mic_offload.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.MPI$(SFX)

uninstall.parallel: clean
	-rm -f $(BINDIR)/pmemd.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.hybrid.MPI$(SFX)
	-rm -f $(BINDIR)/pmemd.hybrid.MPI.mic$(SFX)

depends:
	$(MAKE) -C src/ depends

clean:
	$(MAKE) -C src/ clean

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
	@(if grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type mic_native, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -mic_native." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type mic_offload, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -mic_offload." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(MIC2_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type mic2, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -mic2." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(HYBRID_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type hybrid, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -hybrid." ;\
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

configured_mic: configured
	@(grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic_native!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic_native flag." ;\
		exit 2 ;\
		} ;\
	)

configured_mic2: configured
	echo
	@(grep $(MIC2_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic2!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic2 flag." ;\
		exit 2 ;\
		} ;\
	)

configured_mic_offload: configured 
	echo
	@(grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic_offload!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic_offload flag." ;\
		exit 2 ;\
		} ;\
	)

configured_hybrid: configured 
	echo
	@(grep $(HYBRID_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type hybrid !" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -hybrid flag." ;\
		exit 2 ;\
		} ;\
	)

