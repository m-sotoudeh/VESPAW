# Installing VESPAW

VESPAW uses the same `configure`/`make` build system as upstream CP-PAW.
This page covers the real, tested path; for anything toolchain-specific
not mentioned here, `doc/install.pdf` is still the underlying reference
(the build system itself is unchanged CP-PAW).

## 0. One command (recommended)

```sh
./install.sh                     # default: ~/.vespaw, serial + MPI, with a check
./install.sh --help              # all options
```

What it does:

1. **Toolchain.** `--toolchain auto` (default) uses Homebrew on macOS when
   it is present (`brew install gcc openblas fftw hdf5 open-mpi` for
   anything missing) and a private conda-forge environment (micromamba,
   gfortran 15, OpenBLAS, FFTW, HDF5, Open MPI) everywhere else. No
   administrator rights are needed; nothing outside `--prefix` is changed
   except Homebrew packages. `--toolchain system` (or `--system`) uses
   what is already loaded, for clusters with environment modules.
2. **MPI check.** A two-task test program is compiled and run with the
   chosen MPI before the MPI build. If it fails, only the serial version is
   built and the reason is written to `~/.vespaw/install.log`.
3. **Build.** Writes `parms.auto` for the detected toolchain, runs
   `configure`, `make fast_new`, `make fast_parallel_new` and
   `make scripts`.
4. **Launcher.** Installs `~/.vespaw/bin/vespaw`, which sets up the paths
   and calls the dispatcher; `vespaw doctor` checks the installation
   (executables, libraries, architecture, MPI, setup library, Python).
5. **Check.** Runs the regression cases `pbe_atom` and `si8_gamma`
   (and `si8_gamma` on two MPI tasks) and stops with an error if they fail.

The sections below describe the manual build with a hand-written parameter
file.

## 1. Prerequisites

A Fortran compiler, BLAS/LAPACK, FFTW3, and HDF5. Two toolchains are
known to work end-to-end and have a ready-made parameter file in the
repo root:

| Platform | Compiler/libraries | Parameter file |
|---|---|---|
| macOS, Apple Silicon (arm64) — **default** | native gfortran + OpenBLAS + FFTW3 (Homebrew `/opt/homebrew`) | `parms.in_use` (identical to `parms.arm64`) |
| macOS or Linux, x86_64 | gfortran + OpenBLAS + FFTW3 (Homebrew `/usr/local` or a Linux package manager) | `parms.x86_64` |

On Debian/Ubuntu, `installation/paketliste_debian.txt` lists the system
packages known to cover this.

**Apple Silicon note**: since 2026-09-23 the default `parms.in_use` is
the native arm64 configuration. The x86_64 file (`parms.x86_64`) still
works on an Apple Silicon Mac, but the resulting binary runs under
Rosetta 2 translation and is ~2.2× slower (64-atom Si: 21.6 vs 10.0
min). The arm64-only Davidson divergence that used to block this build
was fixed on 2026-09-23 (see [KNOWN_ISSUES.md](KNOWN_ISSUES.md)).
On x86_64 Linux/macOS, configure with `parms.x86_64`.

## 2. Configure and build

```sh
./configure --with-parmfile=parms.in_use   # Apple Silicon; parms.x86_64 on x86_64
make fast_new                # serial optimized build -> bin/<arch>/vespaw_fast.x
make fast_parallel_new       # optional MPI build (requires mpif90) -> bin/<arch>/vespaw_fast_mpi.x
make scripts                 # installs bin/<arch>/vespaw, the CLI dispatcher (see USAGE.md)
```

`<arch>` is whatever `ARCH=` says inside the parameter file you passed
(`gfortran_openblas_fftw3_arm64` for `parms.in_use`/`parms.arm64`,
`gfortran_openblas_fftw3` for `parms.x86_64`) — each toolchain gets
its own `bin/<arch>/` directory, so you can keep more than one build
around at once.

Put `bin/<arch>/` on your `PATH` (or reference it directly) to run
`vespaw` — see [USAGE.md](USAGE.md).

## 3. Gotchas actually hit building this repo

- **`Permission denied` running `f90pp`** — `config.status` regenerates
  `bin/<arch>/f90pp` without its executable bit. Run
  `chmod +x bin/<arch>/f90pp` after any `./configure` or `./config.status`
  rerun.
- **An edit to `Makefile_bare.in` doesn't take effect** — `./config.status`
  alone only regenerates the top-level `Makefile_bare`/`Makefile`, not the
  per-mode `bin/<arch>/Objects/<mode>/Makefile` copies the build actually
  uses. After editing `Makefile_bare.in`, rerun the full `./configure
  --with-parmfile=...` (not just `./config.status`).
- **`mpirun -np N vespaw ...` looks *slower* than serial** — almost
  always missing `OPENBLAS_NUM_THREADS`. Each MPI rank's own BLAS calls
  (orthogonalization, scalar products) otherwise each spawn a
  full-core-count thread pool independently, oversubscribing the machine
  by roughly `N×`. Set it explicitly:
  ```sh
  OPENBLAS_NUM_THREADS=1 mpirun -np 4 vespaw foo.cntl
  ```
  For a serial run using OpenMP threading instead, set `OMP_NUM_THREADS`
  (and leave `OPENBLAS_NUM_THREADS` unset or matched to it) — see
  `the dev notes`'s Phase 1 profiling notes for the full story of
  how this was found.
- **A fresh build segfaults compiling `vespaw_edft.f90`** (or a couple of
  other specific files) — a pre-existing `f90pp_tmplts.x` bug, already
  worked around for the files known to hit it via a `-notemplates`
  override rule in `Makefile_bare.in`. If a *new* file trips it, add the
  same override rule (search `Makefile_bare.in` for `-notemplates` to see
  the existing pattern) rather than debugging `f90pp_tmplts.x` itself.

## 4. Verifying the build

Once built, run the regression suite (small, fast, checked into the repo
— distinct from the gitignored exploratory `tests/` tree):

```sh
python3 testsuite/run.py
```

See [USAGE.md](USAGE.md) for running your own calculation.
