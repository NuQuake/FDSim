# FDSim3D

Fortran 95 code for finite-difference simulation of seismic wave propagation and earthquake
ground motion in 3D heterogeneous viscoelastic media with a planar free surface.

- Velocity–stress staggered-grid FD scheme, 4th order in space, 2nd order in time
- Discontinuous grid (finer upper part, ratio `REFINE` = any odd number)
- Viscoelastic rheology: generalized Maxwell body ≡ generalized Zener body (arbitrary Q(ω))
- Non-reflecting boundaries, optionally combined with PML
- Excitation: double-couple point sources or vertically incident plane wave
- Serial or MPI (domain decomposition) execution

Authors: Jozef Kristek, Peter Moczo (Comenius University in Bratislava); FDMod3D co-authored by
Martin Gális.

## Repository layout

```
.
├── src/
│   ├── fdmod3d/      FDMod3D    – builds the FD grid model (material cells, effective media, Q)
│   ├── fdmodpar3d/   FDModPar3D – splits the FDMod3D output per MPI process (serial)
│   ├── fdsim3d/      FDSim3D    – FD time-stepping solver
│   └── sourcetf/     SourceTF   – source-time-function generator (serial)
├── examples/
│   └── can4/         sample problem (Can4 basin model, double-couple source, 326 receivers)
├── docs/
│   └── User_Guide_for_FDSim3D.pdf
├── CMakeLists.txt    CMake build
├── Makefile          plain-make build (configuration in make.inc)
├── make.inc
├── CITATION.cff
├── CHANGELOG.md
└── LICENSE           EUPL v1.2
```

## Requirements

- Fortran 95 compiler with a C-style preprocessor (`-cpp` / `-fpp`); tested with GNU Fortran 13
- MPI (MPI-1.2 or newer; OpenMPI and MPICH tested by the authors) for parallel runs
- CMake ≥ 3.18 (optional; plain `make` works as well)

No other libraries are needed.

## Building

### CMake

```bash
cmake -S . -B build                                 # serial, Release
cmake -S . -B build -DFDSIM_ENABLE_MPI=ON           # MPI
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug        # runtime checks, backtraces
cmake --build build -j
# executables: build/bin/{FDMod3D,FDModPar3D,FDSim3D,SourceTF}
```

| Option | Default | Meaning |
|---|---|---|
| `FDSIM_ENABLE_MPI` | `OFF` | build FDMod3D and FDSim3D with MPI (`-DUSE_MPI`) |
| `FDSIM_MPI_MODULE` | `ON` | `USE MPI` (`-DMPI2`); `OFF` → `INCLUDE 'mpif.h'` |
| `FDMOD3D_MOD_FUNC` | `src/fdmod3d/mod_func.F90` | model-definition module compiled into FDMod3D |

Compiler selection: `FC=ifx cmake -S . -B build ...`.

### make

```bash
make                         # all programs, serial, gfortran
make MPI=1                   # FDMod3D + FDSim3D with MPI (mpif90)
make COMPILER=intel MPI=1    # ifort / mpiifort      (COMPILER=ifx → ifx / mpiifx)
make DEBUG=1
make fdsim3d MPI=1           # single program
make clean | distclean
# executables: bin/
```

Further variables (`FC`, `FFLAGS`, `LDFLAGS`, `MPI_MODULE`, `MOD_FUNC`) are documented in `make.inc`.
Objects go to `build-make/<program>-<serial|mpi>[-debug]/`, so serial and MPI builds do not mix;
`bin/` holds the most recently built variant.

As recommended in the user guide, check highly optimised executables against a low-optimisation
(`DEBUG=1` / `Debug`) build.

## Defining a model: `mod_func.F90`

FDMod3D obtains vp, vs, ρ, Qp, Qs from functions in `mod_func.F90`, which has to be written for each
model (User Guide, §4.2). The file in `src/fdmod3d/` defines the Can4 sample model. Keep your own
version next to the model's input files and build against it:

```bash
cmake -S . -B build -DFDMOD3D_MOD_FUNC=$PWD/my_model/mod_func.F90
# or
make fdmod3d MOD_FUNC=$PWD/my_model/mod_func.F90
```

## Running

All programs read and write in the current directory; `HF_FDSim3D` holds the job name `JOBNAME`.

Serial:

1. Prepare `HF_FDSim3D`, `JOBNAME.MD`, interface-surface files, and `mod_func.F90`
2. `FDMod3D` → `JOBNAME.INM`, model files (`*.MO`, `*.Q`, `*.JMH`, `*.JMHF`)
3. Prepare `JOBNAME.IN`, source file (`*.SR`), `SRC_xxx.DAT` (e.g. from `SourceTF`), receiver file
4. `FDSim3D` → `JOBNAME.LOG`, seismograms (`*.asc`), snapshots, `ENERGY.DAT`

Parallel:

```bash
mpirun -np <N1> FDMod3D      # N1 is independent of the FDSim3D process count
FDModPar3D                   # serial; partitions the model for PROCX*PROCY*PROCZ processes
mpirun -np <PROCX*PROCY*PROCZ> FDSim3D
```

Input/output formats are described in `docs/User_Guide_for_FDSim3D.pdf`.

## Example

`examples/can4/` – see [examples/can4/README.md](examples/can4/README.md).

## Citation

Results obtained with the code must reference (see also `CITATION.cff`):

- Kristek, J., Moczo, P., 2014. FDSim3D – The Fortran95 Code for Numerical Simulation of Seismic Wave
  Propagation in 3D Heterogeneous Viscoelastic Media. www.cambridge.org/moczo
- Moczo, P., Kristek, J., Gális, M., 2014. *The Finite-Difference Modelling of Earthquake Motions:
  Waves and Ruptures.* Cambridge University Press.

## License

European Union Public Licence v1.2 – see [LICENSE](LICENSE). Version 1.0 was originally released
in 2014 under EUPL v1.1.
