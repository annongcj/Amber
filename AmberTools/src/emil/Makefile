include ../config.h
include emilfiles

OBJECTS=$(EMIL_CPPSOURCES:.C=.o)

all: install
yes: install
install: $(LIBDIR)/libemil.a

# Only build the static library: the executable
# which used to serve as a test harness for the lib
# hasn't been needed for a while. 
$(LIBDIR)/libemil.a: $(OBJECTS) ../config.h
	@echo "[EMIL]  AR $@"
	$(VB)$(AR) $(LIBDIR)/libemil.a $(OBJECTS)  > /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libemil.a

# All objects depend on all headers.
$(OBJECTS): $(EMIL_CPPHEADERS) ../config.h

.C.o:  
	@echo "[EMIL]  CXX $<"
	$(VB)$(CXX) -c $(CXXOPTFLAGS) $(EMIL_MPIFLAGS) -o $@ $<

.PHONY:clean
clean:
	rm -f *.o *~ \#* 
