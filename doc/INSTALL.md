# Installing VESPAW

## 1. One command (recommended)

```sh
./install.sh                     # default: ~/.vespaw, serial + MPI, with a check
./install.sh --help              # all options
export PATH="$HOME/.vespaw/bin:$PATH"
vespaw doctor                    # verify the installation
```

What `install.sh` does:

1. **Toolchain.** `--toolchain auto` (default) uses Homebrew on macOS when
   it is present (`brew install gcc openblas fftw hdf5 open-mpi` for
   anything missing) and a private conda-forge environment (micromamba,
   gfortran, OpenBLAS, FFTW, HDF5, Open MPI) everywhere else. No
   administrator rights are needed; nothing outside `--prefix` is changed
   except Homebrew packages. `--toolchain system` (or `--system`) uses
   what is already loaded, for clusters with environment modules.
2. **MPI check.** A two-task test program is compiled and run with the
   chosen MPI before the MPI build. If it fails, only the serial version
   is built, and the reason is written to `~/.vespaw/install.log`.
3. **Build.** Writes a parameter file `parms.auto` for the detected
   toolchain, runs `configure`, `make fast_new`, `make fast_parallel_new`
   and `make scripts`.
4. **Launcher.** Installs `~/.vespaw/bin/vespaw`, which sets up the paths
   and calls the dispatcher. `vespaw doctor` checks the executables,
   libraries, architecture, MPI, setup library and Python package.
5. **Check.** Runs the regression cases `pbe_atom` and `si8_gamma` (and
   `si8_gamma` on two MPI tasks) and stops with an error if they fail.

Other options: `--prefix DIR`, `--no-mpi`, `--no-test`, and
`--stage env|build|install|test` to repeat a single stage.

## 2. Manual build

Prerequisites: a Fortran compiler (gfortran 13 to 15 are tested),
BLAS/LAPACK (OpenBLAS), FFTW3 and HDF5 with the Fortran interface; for
the parallel version an MPI library with `mpif90`.

Parameter files for `configure` are in `arch/`:

| Platform | Compiler and libraries | Parameter file |
|---|---|---|
| macOS, Apple Silicon | gfortran, OpenBLAS, FFTW3, HDF5 from Homebrew (`/opt/homebrew`) | `arch/parms.arm64` |
| macOS (Intel) or Linux, x86_64 | gfortran, OpenBLAS, FFTW3, HDF5 (Homebrew `/usr/local` or the system package manager) | `arch/parms.x86_64` |
| macOS, Intel compilers | ifort and MKL | `arch/parms.osx_ifort_mkl` |

Copy the file closest to your system, adjust the library paths, then

```sh
./configure --with-parmfile=arch/parms.x86_64
make fast_new             # serial optimized build  -> bin/<arch>/vespaw_fast.x
make fast_parallel_new    # MPI build (needs mpif90) -> bin/<arch>/vespaw_fast_mpi.x
make scripts              # the vespaw command        -> bin/<arch>/vespaw
```

`<arch>` is the `ARCH=` value of the parameter file, so several builds
can coexist. Put `bin/<arch>/` on your `PATH`, and set
`VESPAW_POT_DIR` to the `potentials/` directory of this repository.

On Debian and Ubuntu the packages `gfortran libopenblas-dev
libfftw3-dev libhdf5-dev libopenmpi-dev openmpi-bin autoconf` are
sufficient; the serial HDF5 Fortran libraries are in
`/usr/lib/x86_64-linux-gnu/hdf5/serial` (add it to `LIBRARY_PATH` and
`LD_LIBRARY_PATH`).

## 3. Known build problems

- **`Permission denied` running `f90pp`.** `config.status` regenerates
  `bin/<arch>/f90pp` without its executable bit. Run
  `chmod +x bin/<arch>/f90pp` after `./configure` or `./config.status`.
- **An edit of `build/Makefile_bare.in` has no effect.** `config.status`
  regenerates only the top-level makefiles, not the copies in
  `bin/<arch>/Objects/<mode>/`. Rerun the full
  `./configure --with-parmfile=...`.
- **`mpirun -np N` is slower than a serial run.** Each MPI task starts its
  own BLAS thread pool, which oversubscribes the machine. Set
  `OPENBLAS_NUM_THREADS=1` (and `OMP_NUM_THREADS=1`) for MPI runs:
  ```sh
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun -np 4 vespaw foo.cntl
  ```
  For a serial run with threads, set both variables to the number of
  cores. Homebrew's OpenBLAS follows `OMP_NUM_THREADS`, the conda-forge
  build follows `OPENBLAS_NUM_THREADS`.
- **The preprocessor crashes on a new source file.** `f90pp_tmplts.x`
  fails on a few files that use no templates; `build/Makefile_bare.in`
  compiles those with `-notemplates`. Add the same rule for the new file.

## 4. Checking the build

```sh
python3 testsuite/run.py pbe_atom si8_gamma   # a few minutes
python3 testsuite/run.py                      # the full suite (22 cases)
VESPAW_REGRESSION_NP=4 python3 testsuite/run.py   # the same with 4 MPI tasks
```

See [USAGE.md](USAGE.md) for running calculations.
