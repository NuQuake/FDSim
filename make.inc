# -----------------------------------------------------------------------------
# Shared build configuration for the plain-Makefile build.
# Included by src/*/Makefile. Override any variable on the command line, e.g.
#
#   make                          # serial, GNU Fortran, optimised
#   make MPI=1                    # MPI build (FDMod3D, FDSim3D)
#   make COMPILER=intel MPI=1     # Intel classic (ifort / mpiifort)
#   make COMPILER=ifx MPI=1       # Intel LLVM (ifx / mpiifx)
#   make DEBUG=1                  # -O0, runtime checks, backtraces
#   make FC=my_f90 FFLAGS="..."   # full manual control
# -----------------------------------------------------------------------------

COMPILER   ?= gnu
MPI        ?= 0
# MPI_MODULE=1 -> "USE MPI"   (-DMPI2)
# MPI_MODULE=0 -> "INCLUDE 'mpif.h'"
MPI_MODULE ?= 1
DEBUG      ?= 0

ifeq ($(COMPILER),gnu)
  SERIAL_FC := gfortran
  MPI_FC    := mpif90
  FPP       := -cpp
  OPTFLAGS  := -O3
  DBGFLAGS  := -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
  # lanczos.F90 uses a cpp line splice ("\") that exceeds 132 columns after preprocessing
  STDFLAGS  := -ffree-line-length-none
  MODFLAG    = -J$(O)
else ifeq ($(COMPILER),intel)
  SERIAL_FC := ifort
  MPI_FC    := mpiifort
  FPP       := -fpp
  OPTFLAGS  := -Ofast
  DBGFLAGS  := -O0 -g -check all -traceback -fpe0
  STDFLAGS  :=
  MODFLAG    = -module $(O)
else ifeq ($(COMPILER),ifx)
  SERIAL_FC := ifx
  MPI_FC    := mpiifx
  FPP       := -fpp
  OPTFLAGS  := -O3
  DBGFLAGS  := -O0 -g -check all -traceback -fpe0
  STDFLAGS  :=
  MODFLAG    = -module $(O)
else
  $(error Unknown COMPILER=$(COMPILER); use gnu, intel or ifx, or set FC/FFLAGS/MODFLAG manually)
endif

ifeq ($(DEBUG),1)
  LEVELFLAGS := $(DBGFLAGS)
  VARIANT_DBG := -debug
else
  LEVELFLAGS := $(OPTFLAGS)
  VARIANT_DBG :=
endif

# Programs that support MPI set USES_MPI := 1 before including this file.
ifeq ($(USES_MPI)$(MPI),11)
  ifeq ($(origin FC),default)
    FC := $(MPI_FC)
  endif
  DEFINES := -DUSE_MPI
  ifeq ($(MPI_MODULE),1)
    DEFINES += -DMPI2
  endif
  VARIANT := mpi$(VARIANT_DBG)
else
  ifeq ($(origin FC),default)
    FC := $(SERIAL_FC)
  endif
  DEFINES :=
  VARIANT := serial$(VARIANT_DBG)
endif

FFLAGS  ?= $(LEVELFLAGS) $(STDFLAGS)
LDFLAGS ?= $(LEVELFLAGS)

# Output locations (relative to src/<program>/)
TOPDIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BINDIR ?= $(TOPDIR)/bin
O      := $(TOPDIR)/build-make/$(PROG)-$(VARIANT)
