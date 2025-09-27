SHELL=/bin/sh
AMBER=Amber24
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA
MIC_TOKEN=mmic
MIC2_TOKEN=DMIC2
MIC_OFFLOAD_TOKEN=DMIC_offload

include $(CONFIG_FILE)

.NOTPARALLEL:configured_cuda configured_parallel configured_serial

install: $(INSTALLTYPE)

pmemd_only:
	cd pmemd && $(MAKE) $(INSTALLTYPE)

cray_serial: serial

cray_parallel: parallel

cray_openmp: openmp

serial: configured_serial
	@echo "Starting installation of ${AMBER} (serial) at `date`".
	cd pmemd && $(MAKE) install
	@echo ""
	@echo "Installation of ${AMBER} (serial) is complete at `date`."
	@echo ""

parallel: configured_parallel
	@echo "Starting installation of ${AMBER} (parallel) at `date`".
	cd ../AmberTools/src/etc && $(MAKE) parallel
	cd pmemd && $(MAKE) parallel
	@echo ""
	@echo "Installation of ${AMBER} (parallel) is complete at `date`."
	@echo ""

cuda_serial: configured_cuda
	@echo "Starting installation of ${AMBER} (cuda serial) at `date`".
	cd pmemd && $(MAKE) cuda_serial

cuda_parallel: configured_cuda configured_parallel
	@echo "Starting installation of ${AMBER} (cuda parallel) at `date`".
	cd ../AmberTools/src/etc && $(MAKE) parallel
	cd pmemd && $(MAKE) cuda_parallel

mic: configured_mic
	@echo "Starting installation of ${AMBER} (mic_native) at `date`".
	cd pmemd && $(MAKE) mic
	@echo ""
	@echo "Installation of ${AMBER} (mic_native) is complete at `date`."
	@echo ""

mic_parallel: configured_mic configured_parallel
	@echo "Starting installation of ${AMBER} (parallel mic_native) at `date`".
	cd ../AmberTools/src/etc && $(MAKE) parallel
	cd pmemd && $(MAKE) mic_parallel
	@echo ""
	@echo "Installation of ${AMBER} (parallel mic_native) is complete at `date`."
	@echo ""

mic2: configured_mic2
	@echo "Starting installation of ${AMBER} (mic2) at `date`".
	cd pmemd && $(MAKE) mic2
	@echo ""
	@echo "Installation of ${AMBER} (mic2) is complete at `date`."
	@echo ""

mic_offload: configured_mic_offload
	@echo "Starting installation of ${AMBER} (mic_offload) at `date`".
	cd pmemd && $(MAKE) mic_offload
	@echo ""
	@echo "Installation of ${AMBER} (mic_offload) is complete at `date`."
	@echo ""

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
	@(if grep $(MIC2_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type mic2, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -mic2." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type mic_offload, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) without -mic_offload." ;\
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
	@(grep $(MIC2_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic2!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic2 flag." ;\
		exit 2 ;\
		} ;\
	)

configured_mic_offload: configured
	@(grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null || \
	{ echo "Error: $(CONFIG_FILE) is not of type mic_offload!" ;\
		echo "  Rerun $(CONFIG_COMMAND) with -mic_offload flag." ;\
		exit 2 ;\
		} ;\
	)

netcdf_uninstall:
	-(cd ../AmberTools/src && $(MAKE) netcdf_uninstall)

clean::
	-(cd pmemd && $(MAKE) clean)

silentclean:
	@$(MAKE) clean 1> /dev/null 2>&1

uninstall: uninstall.serial uninstall.parallel uninstall.cuda uninstall.mic netcdf_uninstall
	-(cd ../test && $(MAKE) clean)
	-/bin/rm -f ../test/numprocs config.h

superclean: uninstall

openmp:
	@echo "No more OpenMP-enabled programs to install"

ompmpi:
	@echo "No more OpenMP/MPI program to install"

uninstall.serial:
	-(cd pmemd && $(MAKE) uninstall.serial)

uninstall.parallel:
	-(cd pmemd && $(MAKE) uninstall.parallel)

uninstall.cuda:
	-(cd pmemd && $(MAKE) uninstall.cuda)

uninstall.mic:
	-(cd pmemd && $(MAKE) uninstall.mic)
