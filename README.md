## Gravothermal SIDM C Solver

This C code solves the gravothermal evolution of self-interacting dark matter halos. It is based on the public Python implementation presented in https://github.com/kboddy/GravothermalSIDM.



## Additional options for initial conditions:
Exponentially truncated  NFW/Hernquist


## Quick Start

Build the executable:

```bash
make
```

Run the default halo:

```bash
make run
```

Run a truncated NFW halo with command-line overrides:

```bash
make run ARGS="--profile tNFW --r-cut 30 --n-shells 800 --t-end 200 --run-name tnfw-test"
```

Run a truncated Hernquist halo:

```bash
make run ARGS="--profile tHern --r-cut 30 --n-shells 800 --t-end 200 --run-name thern-test"
```

Results are written to:

```text
output_dir/run_name/
```

With the default settings this is:

```text
LensingPerturber/8e6-c45-s50-200/
```

Clean the executable:

```bash
make clean
```

## Overview

This directory contains a minimal make-based C implementation of a spherically
symmetric gravothermal evolution solver for self-interacting dark matter
halos. The code evolves a halo in fixed Lagrangian mass coordinates. The shell
radii, density, pressure, and derived thermal quantities are updated through
alternating heat-conduction and hydrostatic-adjustment steps.

The implementation is intentionally plain:

- `gravothermal_sidm.c`: main implementation and executable entry point
- `gravothermal_sidm.h`: configuration structs and public function declarations
- `Makefile`: single-command build and run targets

No CMake project or build directory is required. Running `make` creates one
executable named `gravothermal_sidm` in this directory.

## Requirements

Required:

- A C11 compiler, usually `gcc` or `clang`
- `make`
- The standard math library, linked with `-lm`

Optional:

- OpenMP support for parallel shell-local loops

The default Makefile enables OpenMP:

```bash
make USE_OPENMP=1
```

If your compiler does not support OpenMP, build without it:

```bash
make clean
make USE_OPENMP=0
```

## Build And Run

Compile only:

```bash
make
```

Compile and run:

```bash
make run
```

Pass runtime arguments through `ARGS`:

```bash
make run ARGS="--run-name myrun --output-dir results --n-shells 1000 --t-end 500"
```

Run the executable directly:

```bash
./gravothermal_sidm --run-name myrun --output-dir results
```

Show all command-line options:

```bash
./gravothermal_sidm --help
```

## Default Parameters

The default run is set in:

```text
gravothermal_sidm.c
function: gsidm_set_default_config()
```

Edit that function if you want to change the default run without passing
command-line arguments. Most run settings have command-line overrides; a few
model constants such as `a`, `b`, `C`, and `gamma` are currently edited in the
source.

The shipped defaults mirror the original `run8e6-c45-s50-200.py` runner:

```text
run_name = 8e6-c45-s50-200
output_dir = LensingPerturber
profile = NFW
r_s = 0.09265 kpc
rho_s = 0.2808 Msun/pc^3
sigma/m = 50 cm^2/g
n_shells = 200
r_min = 0.005
r_max = 60
t_end = 2000
```

## Command-Line Parameters

General run control:

| Option | Meaning | Unit |
|---|---|---|
| `--run-name NAME` | Name of this run; also the output subdirectory | text |
| `--output-dir DIR` | Parent output directory | path |
| `--resume latest` | Resume from the latest snapshot in `output_dir/run_name` | text |
| `--resume FILE` | Resume from a specific snapshot file | path or filename |
| `--help` | Print command-line help | none |

Initial halo parameters:

| Option | Meaning | Unit |
|---|---|---|
| `--profile NAME` | Initial profile: `NFW`, `Hernquist`, `Isothermal`, `tNFW`, or `tHern` | text |
| `--r-s VALUE` | Scale radius `r_s` | kpc |
| `--rho-s VALUE` | Scale density `rho_s` | Msun/pc^3 |
| `--r-cut VALUE` | Exponential cutoff radius for `tNFW` and `tHern` | dimensionless, in units of `r_s` |
| `--sigma-m-with-units VALUE` | Self-interaction cross section per unit mass | cm^2/g |
| `--w-units VALUE` | Velocity scale used by velocity-dependent cross-section models | km/s |
| `--n-shells VALUE` | Number of radial shells | integer |
| `--r-min VALUE` | Minimum initial radius | dimensionless, `r/r_s` |
| `--r-max VALUE` | Maximum initial radius | dimensionless, `r/r_s` |
| `--p-rmax-factor VALUE` | Pressure integration upper limit is `p_rmax_factor*r_max` | dimensionless |

Numerical evolution parameters:

| Option | Meaning | Unit |
|---|---|---|
| `--t-end VALUE` | Stop when dimensionless time reaches this value | dimensionless |
| `--rho-factor-end VALUE` | Stop when central density reaches this factor times its initial value | dimensionless |
| `--kn-end VALUE` | Stop when the minimum Knudsen number drops below this value; `0` disables it | dimensionless |
| `--t-epsilon VALUE` | Timestep safety factor | dimensionless |
| `--r-epsilon VALUE` | Hydrostatic-adjustment tolerance for relative radius change | dimensionless |
| `--n-adjustment-fixed VALUE` | Force this many hydrostatic-adjustment iterations per timestep; `-1` disables | integer |
| `--n-adjustment-max VALUE` | Maximum hydrostatic-adjustment iterations per timestep; `-1` disables | integer |
| `--flag-timestep-use-relaxation 0|1` | Use relaxation-time timestep criterion | boolean |
| `--flag-timestep-use-energy 0|1` | Use energy-change timestep criterion | boolean |
| `--flag-hydrostatic-initial 0|1` | Apply hydrostatic adjustment to the initial condition | boolean |

Snapshot cadence:

| Option | Meaning | Unit |
|---|---|---|
| `--save-frequency-rate VALUE` | Save based on `t/t_relax`; use `none` to disable | dimensionless |
| `--save-frequency-timing VALUE` | Save when relative time increase exceeds this value; use `none` to disable | dimensionless |
| `--save-frequency-density VALUE` | Save when central density changes by this fraction; use `none` to disable | dimensionless |

OpenMP options:

| Option | Meaning | Unit |
|---|---|---|
| `--use-openmp 0|1` | Enable OpenMP regions at runtime if compiled with OpenMP | boolean |
| `--threads VALUE` | Number of OpenMP threads | integer |

Scattering model options:

| Option | Meaning |
|---|---|
| `--lmfp-model NAME` | Long-mean-free-path scattering model |
| `--smfp-model NAME` | Short-mean-free-path scattering model |

Currently supported scattering model strings include:

```text
constant
powerlaw_INDEX
YukawaBornViscosityApprox_K3_order1
YukawaBornViscosityApprox_K5_order1
YukawaBornViscosityApprox_K7_order1
YukawaBornViscosityApprox_K9_order1
```

## Initial Density Profiles

All profiles are expressed in the dimensionless radius:

```text
x = r / r_s
```

The density arrays stored by the code are dimensionless:

```text
rho_dimensionless = rho_physical / rho_s
```

Supported profiles:

```text
NFW:
rho(x) = 1 / [x (1 + x)^2]

Hernquist:
rho(x) = 1 / [x (1 + x)^3]

Isothermal:
rho(x) = 1 / (1 + x^2)

tNFW:
rho(x) = exp(-x / r_cut) / [x (1 + x)^2]

tHern:
rho(x) = exp(-x / r_cut) / [x (1 + x)^3]
```

For the truncated profiles, `r_cut` is dimensionless. If the physical cutoff
radius is `R_cut`, pass:

```text
--r-cut R_cut/r_s
```

The code initializes pressure from hydrostatic equilibrium. For `NFW`, `tNFW`,
and `tHern`, the enclosed mass and pressure integrals are evaluated
numerically during initialization. For `Hernquist` and `Isothermal`, analytic
pressure expressions are used.

## Internal Units

The code evolves dimensionless variables. The conversion factors are written to
`halo_init.txt` in each output directory.

The fundamental input scales are:

```text
scale_r   = r_s
scale_rho = rho_s
scale_m   = 4 pi rho_s r_s^3
```

In SI units, the code writes:

```text
scale_r_m
scale_rho_kg_per_m3
scale_m_kg
scale_u_m2_per_s2
scale_p_pa
scale_v_m_per_s
scale_t_s
scale_L_w
scale_sigma_m_m2_per_kg
```

### Converting Snapshot Arrays

Let `r`, `m`, `rho`, `p`, and `Mcy` be the arrays read from a binary snapshot.

Radius:

```text
R_m   = r * scale_r_m
R_kpc = r * r_s_kpc
```

Density:

```text
rho_kg_m3      = rho * scale_rho_kg_per_m3
rho_Msun_pc3   = rho * rho_s_msun_per_pc3
```

Spherical enclosed mass:

```text
M_kg   = m * scale_m_kg
M_Msun = m * scale_m_kg / M_sun
```

Pressure:

```text
P_pa = p * scale_p_pa
```

Projected cylindrical mass:

```text
Mcy_kg   = Mcy * scale_m_kg
Mcy_Msun = Mcy * scale_m_kg / M_sun
```

One-dimensional velocity dispersion:

```text
sigma_1D_dimensionless = sqrt(p / rho)
sigma_1D_m_s           = sqrt(p / rho) * scale_v_m_per_s
sigma_1D_km_s          = sqrt(p / rho) * scale_v_m_per_s / 1000
```

This follows the code definition:

```text
v = sqrt(p/rho)
u = 3 p / (2 rho)
```

Time:

```text
t_s   = t * scale_t_s
t_Gyr = t * scale_t_s / (3600 * 24 * 365.25 * 1e9)
```

## Output Files

Each run writes files under:

```text
output_dir/run_name/
```

Main files:

| File | Meaning |
|---|---|
| `halo_init.txt` | Run parameters and unit-conversion factors |
| `snapshot_index.csv` | One row per saved snapshot |
| `time_*.bin` | Binary state snapshots |

### `halo_init.txt`

This is a text file with `key=value` lines. It records the resolved runtime
configuration and all unit conversion factors. Important entries include:

```text
r_s_kpc
rho_s_msun_per_pc3
r_cut
sigma_m_with_units_cm2_per_g
scale_r_m
scale_rho_kg_per_m3
scale_m_kg
scale_p_pa
scale_v_m_per_s
scale_t_s
```

### `snapshot_index.csv`

Columns:

```text
save_id,prefix,t,t_before,rho_center,n_adjustment,n_conduction,filename
```

Use `save_id` or `t` to order snapshots. Do not rely on lexicographic filename
ordering, because scientific notation and integer strings do not sort in time
order.

### Binary Snapshot Format

Snapshots are written as native C binary data. On normal x86 Linux clusters
this is little-endian. The current format is `version = 2`.

Header and scalar metadata:

| Order | Type | Name |
|---|---|---|
| 1 | `uint64` | magic |
| 2 | `uint64` | version |
| 3 | `uint64` | n_shells |
| 4 | `double` | rho_center |
| 5 | `int64` | n_adjustment |
| 6 | `int64` | n_conduction |
| 7 | `uint64` | n_save |
| 8 | `double` | t_epsilon |
| 9 | `double` | r_epsilon |
| 10 | `double` | t |
| 11 | `double` | t_before |

Then come `n_shells` values for each array, all as `double`:

```text
m
r
rho
p
Mcy
```

Version 1 snapshots from older code contain only:

```text
m
r
rho
p
```

The current code can resume from both version 1 and version 2 snapshots.

## Reading Output In Mathematica

This reader returns `t`, `m`, `r`, `rho`, `p`, and `Mcy`. It also computes the
dimensionless one-dimensional velocity dispersion.

```wolfram
ReadGSIDMSnapshot[file_] :=
 Module[{s, magic, version, n, rhoCenter, nAdjustment, nConduction,
   nSave, tEpsilon, rEpsilon, t, tBefore, m, r, rho, p, mcy},

  s = OpenRead[file, BinaryFormat -> True];

  magic = BinaryRead[s, "UnsignedInteger64", ByteOrdering -> -1];
  version = BinaryRead[s, "UnsignedInteger64", ByteOrdering -> -1];
  n = BinaryRead[s, "UnsignedInteger64", ByteOrdering -> -1];

  rhoCenter = BinaryRead[s, "Real64", ByteOrdering -> -1];
  nAdjustment = BinaryRead[s, "Integer64", ByteOrdering -> -1];
  nConduction = BinaryRead[s, "Integer64", ByteOrdering -> -1];
  nSave = BinaryRead[s, "UnsignedInteger64", ByteOrdering -> -1];
  tEpsilon = BinaryRead[s, "Real64", ByteOrdering -> -1];
  rEpsilon = BinaryRead[s, "Real64", ByteOrdering -> -1];
  t = BinaryRead[s, "Real64", ByteOrdering -> -1];
  tBefore = BinaryRead[s, "Real64", ByteOrdering -> -1];

  m = BinaryReadList[s, "Real64", n, ByteOrdering -> -1];
  r = BinaryReadList[s, "Real64", n, ByteOrdering -> -1];
  rho = BinaryReadList[s, "Real64", n, ByteOrdering -> -1];
  p = BinaryReadList[s, "Real64", n, ByteOrdering -> -1];

  mcy =
   If[version >= 2,
    BinaryReadList[s, "Real64", n, ByteOrdering -> -1],
    Missing["NotAvailable"]
   ];

  Close[s];

  <|
   "file" -> file,
   "version" -> version,
   "n" -> n,
   "rho_center" -> rhoCenter,
   "n_adjustment" -> nAdjustment,
   "n_conduction" -> nConduction,
   "n_save" -> nSave,
   "t_epsilon" -> tEpsilon,
   "r_epsilon" -> rEpsilon,
   "t" -> t,
   "t_before" -> tBefore,
   "m" -> m,
   "r" -> r,
   "rho" -> rho,
   "p" -> p,
   "Mcy" -> mcy,
   "sigma1D" -> Sqrt[p/rho]
   |>
  ]
```

Read all snapshots in a run directory:

```wolfram
ReadAllGSIDM[dir_] :=
 Module[{files, data},
  files = FileNames["time_*.bin", dir];
  data = ReadGSIDMSnapshot /@ files;
  SortBy[data, #["t"] &]
  ]
```

Example:

```wolfram
data = ReadAllGSIDM["LensingPerturber/tnfw-test"];

times = data[[All, "t"]];
rFirst = data[[1, "r"]];
rhoFirst = data[[1, "rho"]];
mcyFirst = data[[1, "Mcy"]];
sigmaFirst = data[[1, "sigma1D"]];
```

To convert `sigma1D` to km/s, multiply by `scale_v_m_per_s/1000` from
`halo_init.txt`.

## Reading Output In Python

This reader assumes little-endian binary files, which is the standard on common
x86 Linux clusters.

```python
from pathlib import Path
import struct
import numpy as np


def read_gsidm_snapshot(path):
    path = Path(path)
    with path.open("rb") as f:
        magic, version, n = struct.unpack("<QQQ", f.read(24))

        rho_center = struct.unpack("<d", f.read(8))[0]
        n_adjustment = struct.unpack("<q", f.read(8))[0]
        n_conduction = struct.unpack("<q", f.read(8))[0]
        n_save = struct.unpack("<Q", f.read(8))[0]
        t_epsilon = struct.unpack("<d", f.read(8))[0]
        r_epsilon = struct.unpack("<d", f.read(8))[0]
        t = struct.unpack("<d", f.read(8))[0]
        t_before = struct.unpack("<d", f.read(8))[0]

        m = np.fromfile(f, dtype="<f8", count=n)
        r = np.fromfile(f, dtype="<f8", count=n)
        rho = np.fromfile(f, dtype="<f8", count=n)
        p = np.fromfile(f, dtype="<f8", count=n)
        Mcy = np.fromfile(f, dtype="<f8", count=n) if version >= 2 else None

    return {
        "file": str(path),
        "version": version,
        "n": n,
        "rho_center": rho_center,
        "n_adjustment": n_adjustment,
        "n_conduction": n_conduction,
        "n_save": n_save,
        "t_epsilon": t_epsilon,
        "r_epsilon": r_epsilon,
        "t": t,
        "t_before": t_before,
        "m": m,
        "r": r,
        "rho": rho,
        "p": p,
        "Mcy": Mcy,
        "sigma1D": np.sqrt(p / rho),
    }


def read_all_gsidm(run_dir):
    files = sorted(Path(run_dir).glob("time_*.bin"))
    data = [read_gsidm_snapshot(file) for file in files]
    return sorted(data, key=lambda item: item["t"])
```

Example:

```python
data = read_all_gsidm("LensingPerturber/tnfw-test")
first = data[0]

r = first["r"]
rho = first["rho"]
Mcy = first["Mcy"]
sigma1D = first["sigma1D"]
```

## Cylindrical Mass `Mcy`

`Mcy[i]` is the projected cylindrical mass enclosed within projected radius:

```text
R = r[i]
```

The continuous definition is:

```text
M_cyl(<R) = 2 pi int_0^R S dS int_-inf^inf rho(sqrt(S^2 + z^2)) dz
```

The code approximates the halo as piecewise constant-density spherical shells.
For shell `j`:

```text
r_{j-1} < r < r_j
rho = rho_j
```

The shell contribution to the dimensionless cylindrical mass is:

```text
Delta Mcy_j(<R) = rho_j/3 * [F(r_j, R) - F(r_{j-1}, R)]
```

where:

```text
F(a, R) = a^3,                         R >= a
F(a, R) = a^3 - (a^2 - R^2)^(3/2),      R < a
```

The total stored value is:

```text
Mcy(<R_i) = sum_j Delta Mcy_j(<R_i)
```

The factor `4 pi` is not written explicitly because the code's dimensionless
mass unit already divides physical mass by:

```text
scale_m = 4 pi rho_s r_s^3
```

Therefore:

```text
M_cyl_physical(<R_i) = Mcy[i] * scale_m_kg
```

## Resume From Snapshot

Resume the latest snapshot for the same `output_dir/run_name`:

```bash
./gravothermal_sidm --resume latest
```

Resume and extend the final time:

```bash
./gravothermal_sidm --resume latest --t-end 3000
```

Resume a specific snapshot:

```bash
./gravothermal_sidm --resume time_25_t_8d150584418em01.bin
```

Through `make`:

```bash
make run ARGS="--resume latest --t-end 3000"
```

Resume behavior:

- The code reads `output_dir/run_name/halo_init.txt`.
- It restores the saved halo configuration.
- It loads the selected binary snapshot.
- It recomputes derived quantities and continues evolving.

Runtime controls such as `--t-end`, save-frequency options, `--use-openmp`, and
`--threads` can still be changed when resuming.

## OpenMP And Performance Notes

The Makefile enables OpenMP by default:

```bash
make USE_OPENMP=1
```

Disable OpenMP:

```bash
make clean
make USE_OPENMP=0
```

Set the shell-count threshold for enabling OpenMP loops:

```bash
make OMP_MIN_SHELLS=4096
```

Runtime threading:

```bash
./gravothermal_sidm --use-openmp 1 --threads 4
```

For small runs such as `n_shells=200`, single-threaded execution is often
fastest. For larger runs, try `--threads 2` or `--threads 4` and benchmark on
your machine. Running many independent halos in parallel is usually more
efficient than using many threads for one small halo.

## Practical Examples

Run a higher-resolution NFW halo:

```bash
make run ARGS="--profile NFW --n-shells 1000 --t-end 500 --run-name nfw-n1000"
```

Run a truncated NFW halo with cutoff at `30 r_s`:

```bash
make run ARGS="--profile tNFW --r-cut 30 --n-shells 1000 --run-name tnfw-rcut30"
```

Run a truncated Hernquist halo and save by density change:

```bash
make run ARGS="--profile tHern --r-cut 30 --save-frequency-density 0.05 --run-name thern-rcut30"
```

Disable all intermediate snapshots except initial and final:

```bash
make run ARGS="--save-frequency-rate none --save-frequency-timing none --save-frequency-density none"
```

Save more frequently by relative time interval:

```bash
make run ARGS="--save-frequency-density none --save-frequency-timing 0.02"
```

## Notes And Caveats

- The radial grid is initialized uniformly in `log(r/r_s)`, but the evolution
  uses fixed Lagrangian mass coordinates. Shell radii move during hydrostatic
  adjustment.
- During deep core collapse, the radial grid can become highly non-uniform.
  Increase `--n-shells` and perform convergence tests for production runs.
- Snapshot files are compact binary files, not portable text files. Use the
  readers above or your own reader matching the documented binary layout.
- The code writes binary data in the machine's native C representation. On
  standard x86 Linux systems this is little-endian IEEE-754 double precision.
