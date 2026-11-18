# Can4 – sample problem

Model Can4 (User Guide, §4.2, Figs. 3–4): a sedimentary basin with three layers over bedrock,
basin edge with linearly varying thickness; double-couple point source at 3 km depth;
326 receivers.

| File | Read by | Content |
|---|---|---|
| `HF_FDSim3D` | all | job name (`CAN4`) |
| `CAN4.MD` | FDMod3D | grid spacing (H = 190 m, REFINE = 19), model boundaries, attenuation frequency band |
| `INTERFACE.DAT` | FDMod3D | material-interface surface |
| `CAN4.IN` | FDSim3D | control data (TW = 50 s, DT = 0.001 s, 8×8×1 MPI processes, PML, snapshots) |
| `CAN45.SR` | FDSim3D | source definition |
| `SRC_001.DAT` | FDSim3D | source time function |
| `RECEIVERS.DAT` | FDSim3D | receiver names and coordinates |

The model itself (vp, vs, ρ, Q) is defined in `src/fdmod3d/mod_func.F90` (the default `MOD_FUNC`).

## Run

```bash
make MPI=1                      # or: cmake -S . -B build -DFDSIM_ENABLE_MPI=ON && cmake --build build
mkdir -p runs/can4 && cp examples/can4/* runs/can4/ && cd runs/can4
mpirun -np 16 ../../bin/FDMod3D
../../bin/FDModPar3D
mpirun -np 64 ../../bin/FDSim3D   # PROCX*PROCY*PROCZ = 8*8*1 from CAN4.IN
```

Resource notes: the finer grid is 1521 × 1521 × 39 cells (h = 10 m). A serial FDMod3D run allocates
~17 GB in `prepare_int` alone, so the example is intended for a cluster. In a serial run of FDSim3D,
`PROCX/PROCY/PROCZ` are ignored and FDModPar3D is not needed.
