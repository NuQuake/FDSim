# FDModPar3D - partitions FDMod3D output for an MPI run of FDSim3D (serial program)
PROG     := FDModPar3D
USES_MPI := 0
include ../../make.inc

SRCS := FDModPar3D.f90
OBJS := $(SRCS:%.f90=$(O)/%.o)

.PHONY: all clean
all: $(O)/$(PROG)
	@mkdir -p $(BINDIR)
	cp -p $< $(BINDIR)/$(PROG)

$(O)/$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

$(O)/%.o: %.f90 | $(O)
	$(FC) $(FFLAGS) $(MODFLAG) -I$(O) -c $< -o $@

$(O):
	mkdir -p $@

clean:
	rm -rf $(O)
