# Changelog

## Unreleased
- Repository restructured for Git hosting: `src/`, `examples/`, `docs/`
- CMake build and unified plain-make build (`make.inc`: serial/MPI, GNU/Intel, debug)
- `mod_func.F90` selectable at build time (`FDMOD3D_MOD_FUNC` / `MOD_FUNC`)
- GitHub Actions compile check (GNU Fortran, serial and OpenMPI)
- Licence upgraded from EUPL v1.1 to EUPL v1.2
- Numerical Recipes code replaced by independent implementations:
  - `fdsim3d/bcucof.F90`: bicubic interpolation coefficients (Hermite form)
  - `fdsim3d/open_sr_file.F90`: clamped cubic spline for source time functions
    (replaces `spline`, `tridag`, `splint`, `locate`)
  - `fdmod3d/nrtype.F90`: reduced to the kinds actually used (SP, DP, WP)
  - `fdmod3d/nrutil.F90`: removed (unused)
  FDMod3D output is bit-identical to version 1.0; FDSim3D seismograms agree to
  within the last printed digit (single-precision rounding).
- Licence headers added to `qylsm.F90`, `write_grd.F90`, `FDModPar3D.f90`, `fcoolr.f90`,
  `mod_prec.f90`, `sourcetf.f90`
- Line endings normalised (CRLF → LF)

## 1.0 – 2014-03-19
- Initial version (Jozef Kristek, Peter Moczo, <fdsim@nuquake.eu>)
