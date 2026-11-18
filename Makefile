# Top-level Makefile: builds all programs into ./bin
# Options (see make.inc): COMPILER=gnu|intel|ifx  MPI=0|1  MPI_MODULE=0|1  DEBUG=0|1  MOD_FUNC=<file>
#
#   make              # all programs, serial
#   make MPI=1        # FDMod3D and FDSim3D with MPI; FDModPar3D, SourceTF stay serial
#   make fdsim3d MPI=1
#   make clean        # remove objects;  make distclean  also removes ./bin

PROGRAMS := fdmod3d fdmodpar3d fdsim3d sourcetf

.PHONY: all $(PROGRAMS) clean distclean
all: $(PROGRAMS)

$(PROGRAMS):
	$(MAKE) -C src/$@

clean:
	rm -rf build-make

distclean: clean
	rm -rf bin
