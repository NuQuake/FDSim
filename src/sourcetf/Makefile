# SourceTF - generates source time functions (serial program)
PROG     := SourceTF
USES_MPI := 0
include ../../make.inc

SRCS := fcoolr.f90 mod_prec.f90 sourcetf.f90
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

include deps.mk
