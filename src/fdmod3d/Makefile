# FDMod3D - grid-model preparation (serial or MPI)
#
# mod_func.F90 defines the medium (vp, vs, rho, Qp, Qs as functions of x,y,z)
# and must be adapted to each model. Point MOD_FUNC to your own file:
#   make MOD_FUNC=/path/to/my_model/mod_func.F90
PROG     := FDMod3D
USES_MPI := 1
include ../../make.inc

MOD_FUNC ?= mod_func.F90

SRCS := $(sort $(filter-out mod_func.F90,$(wildcard *.F90)))
OBJS := $(SRCS:%.F90=$(O)/%.o) $(O)/mod_func.o

.PHONY: all clean
all: $(O)/$(PROG)
	@mkdir -p $(BINDIR)
	cp -p $< $(BINDIR)/$(PROG)

$(O)/$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

$(O)/mod_func.o: $(MOD_FUNC) | $(O)
	$(FC) $(FPP) $(DEFINES) $(FFLAGS) $(MODFLAG) -I$(O) -c $< -o $@

$(O)/%.o: %.F90 | $(O)
	$(FC) $(FPP) $(DEFINES) $(FFLAGS) $(MODFLAG) -I$(O) -c $< -o $@

$(O):
	mkdir -p $@

clean:
	rm -rf $(O)

include deps.mk
