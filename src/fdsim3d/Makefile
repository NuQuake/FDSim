# FDSim3D - finite-difference wave propagation solver (serial or MPI)
PROG     := FDSim3D
USES_MPI := 1
include ../../make.inc

SRCS := $(sort $(wildcard *.F90))
OBJS := $(SRCS:%.F90=$(O)/%.o)

.PHONY: all clean
all: $(O)/$(PROG)
	@mkdir -p $(BINDIR)
	cp -p $< $(BINDIR)/$(PROG)

$(O)/$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

$(O)/%.o: %.F90 | $(O)
	$(FC) $(FPP) $(DEFINES) $(FFLAGS) $(MODFLAG) -I$(O) -c $< -o $@

$(O):
	mkdir -p $@

clean:
	rm -rf $(O)

include deps.mk
