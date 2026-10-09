#!/usr/bin/env bash
# =============================================================================
#  VESPAW one-step installer
#
#    ./install.sh                 build serial + MPI, install, run a check
#    ./install.sh --help          all options
#
#  The installer finds or installs everything VESPAW needs (Fortran
#  compiler, BLAS/LAPACK, FFTW, HDF5, MPI, Python) without administrator
#  rights:
#    macOS with Homebrew   -> Homebrew packages (brew install if missing)
#    otherwise             -> a private conda-forge environment (micromamba)
#    --system              -> what is already loaded (clusters: `module load`)
#  MPI is checked with a two-task test before it is used; if it does not
#  work, the serial version is built and the problem is reported.
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${HOME}/.vespaw"
WITH_MPI=1
RUN_TESTS=1
TOOLCHAIN="auto"
STAGE="all"

usage() {
  cat <<EOF
usage: ./install.sh [options]

  --prefix DIR        install location (default: ~/.vespaw)
  --toolchain NAME    auto (default) | brew | conda | system
                        brew   : Homebrew (macOS); missing packages are installed
                        conda  : private conda-forge environment in PREFIX/env
                        system : compilers/libraries already in PATH (clusters:
                                 load gfortran, MPI, OpenBLAS, FFTW, HDF5 first)
  --system            same as --toolchain system
  --no-mpi            build the serial version only
  --no-test           skip the post-install check
  --stage STAGE       run only one stage: env | build | install | test
  -h, --help          this help

After installation, put PREFIX/bin on your PATH; then
  vespaw FILE.cntl        run a calculation
  vespaw doctor           check the installation
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2 ;;
    --prefix=*) PREFIX="${1#*=}"; shift ;;
    --toolchain) TOOLCHAIN="$2"; shift 2 ;;
    --toolchain=*) TOOLCHAIN="${1#*=}"; shift ;;
    --system) TOOLCHAIN="system"; shift ;;
    --no-mpi) WITH_MPI=0; shift ;;
    --no-test) RUN_TESTS=0; shift ;;
    --stage) STAGE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1"; usage; exit 1 ;;
  esac
done
mkdir -p "$PREFIX/bin"
PREFIX="$(cd "$PREFIX" && pwd)"
ENV="$PREFIX/env"
LOG="$PREFIX/install.log"
: >> "$LOG"

say()  { printf '\033[1m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33mWARNING: %s\033[0m\n' "$*"; }
fail() { printf '\033[31mERROR: %s\033[0m\n' "$*" >&2; echo "see $LOG" >&2; exit 1; }
want() { [ "$STAGE" = "all" ] || [ "$STAGE" = "$1" ]; }

# -----------------------------------------------------------------------------
#  platform and toolchain choice
# -----------------------------------------------------------------------------
OS="$(uname -s)"; MACH="$(uname -m)"
case "$OS-$MACH" in
  Darwin-arm64)   PLAT=osx-arm64 ;;
  Darwin-x86_64)  PLAT=osx-64 ;;
  Linux-x86_64)   PLAT=linux-64 ;;
  Linux-aarch64)  PLAT=linux-aarch64 ;;
  *) fail "unsupported platform $OS-$MACH" ;;
esac
if [ "$OS" = "Darwin" ]; then
  xcrun --show-sdk-path >/dev/null 2>&1 \
    || fail "macOS needs the Xcode command line tools: run  xcode-select --install"
  export SDKROOT="$(xcrun --show-sdk-path)"
fi
BREW=""
if [ "$OS" = "Darwin" ]; then
  for b in /opt/homebrew/bin/brew /usr/local/bin/brew; do
    # on Apple Silicon only the native (arm64) Homebrew is usable
    if [ -x "$b" ] && { [ "$MACH" != "arm64" ] || [ "$b" = /opt/homebrew/bin/brew ]; }; then
      BREW="$b"; break
    fi
  done
fi
if [ "$TOOLCHAIN" = "auto" ]; then
  if [ -n "$BREW" ]; then TOOLCHAIN=brew; else TOOLCHAIN=conda; fi
fi
say "VESPAW installer: $PLAT, toolchain $TOOLCHAIN, prefix $PREFIX"

# -----------------------------------------------------------------------------
#  1. toolchain
# -----------------------------------------------------------------------------
MAMBA="$PREFIX/bin/micromamba"
get_micromamba() {
  if [ ! -x "$MAMBA" ]; then
    say "downloading micromamba"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${PLAT}/latest" \
      | tar -xj -C "$PREFIX" bin/micromamba || fail "micromamba download failed"
  fi
}
# gfortran 15: VESPAW is tested with it (the build's template preprocessor
# was first tested with 16 in 2026-10 and needed a larger stack there).
CONDA_PKGS="python=3.12 numpy pyyaml h5py make gfortran=15 openblas fftw hdf5"
PY_PKGS="python=3.12 numpy pyyaml h5py"

if want env; then
  case "$TOOLCHAIN" in
    conda)
      get_micromamba
      PKGS="$CONDA_PKGS"; [ "$WITH_MPI" -eq 1 ] && PKGS="$PKGS openmpi"
      say "installing the toolchain from conda-forge into $ENV (a few minutes)"
      if [ -x "$ENV/bin/python3" ]; then
        "$MAMBA" install -y -q -r "$PREFIX/mamba" -p "$ENV" -c conda-forge $PKGS >>"$LOG" 2>&1
      else
        "$MAMBA" create -y -q -r "$PREFIX/mamba" -p "$ENV" -c conda-forge $PKGS >>"$LOG" 2>&1
      fi || fail "toolchain installation failed"
      ;;
    brew)
      [ -n "$BREW" ] || fail "Homebrew not found (install it, or use --toolchain conda)"
      FORMULAE="gcc openblas fftw hdf5"; [ "$WITH_MPI" -eq 1 ] && FORMULAE="$FORMULAE open-mpi"
      MISSING=""
      for f in $FORMULAE; do "$BREW" list --versions "$f" >/dev/null 2>&1 || MISSING="$MISSING $f"; done
      if [ -n "$MISSING" ]; then
        say "installing with Homebrew:$MISSING"
        "$BREW" install $MISSING >>"$LOG" 2>&1 || fail "brew install$MISSING failed"
      else
        say "Homebrew packages present: $FORMULAE"
      fi
      # Python for the tools (YAML input, Python API, regression suite)
      get_micromamba
      if [ ! -x "$ENV/bin/python3" ]; then
        say "installing Python for the VESPAW tools into $ENV"
        "$MAMBA" create -y -q -r "$PREFIX/mamba" -p "$ENV" -c conda-forge $PY_PKGS >>"$LOG" 2>&1 \
          || fail "Python environment failed"
      fi
      ;;
    system) ;;
    *) fail "unknown toolchain $TOOLCHAIN" ;;
  esac
fi

# -----------------------------------------------------------------------------
#  2. locate compilers and libraries
# -----------------------------------------------------------------------------
PY="python3"; [ -x "$ENV/bin/python3" ] && PY="$ENV/bin/python3"
case "$TOOLCHAIN" in
  conda)
    TOOL="$ENV"; FC="$ENV/bin/gfortran"; MPIFC="$ENV/bin/mpifort"; MPIRUN="$ENV/bin/mpirun"
    INCDIRS="-I$ENV/include"; LIBDIRS="-L$ENV/lib -Wl,-rpath,$ENV/lib"
    HDF5LIBS="-lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5"
    FFTH="$ENV/include/fftw3.f03"; MPIH="$ENV/include/mpif.h"
    ;;
  brew)
    BP="$("$BREW" --prefix)"
    TOOL="$BP"; FC="$BP/bin/gfortran"; MPIFC="$BP/bin/mpif90"; MPIRUN="$BP/bin/mpirun"
    OB="$("$BREW" --prefix openblas)"; FW="$("$BREW" --prefix fftw)"; H5="$("$BREW" --prefix hdf5)"
    INCDIRS="-I$FW/include -I$BP/include -I$H5/include"
    LIBDIRS="-L$OB/lib -L$FW/lib -L$H5/lib -L$BP/lib -Wl,-rpath,$OB/lib -Wl,-rpath,$H5/lib"
    HDF5LIBS="-lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5"
    FFTH="$FW/include/fftw3.f03"; MPIH="$BP/include/mpif.h"
    ;;
  system)
    FC="$(command -v gfortran || true)"; MPIFC="$(command -v mpifort || command -v mpif90 || true)"
    MPIRUN="$(command -v mpirun || command -v mpiexec || true)"
    TOOL="$(cd "$(dirname "${FC:-/usr/bin/gfortran}")/.." && pwd)"
    LIBDIRS=""; INCDIRS=""
    for d in ${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH//:/ }} ${LIBRARY_PATH:+${LIBRARY_PATH//:/ }} /usr/lib64 /usr/lib /usr/lib/x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/hdf5/serial; do
      [ -d "$d" ] && LIBDIRS="$LIBDIRS -L$d"
    done
    for d in ${CPATH:+${CPATH//:/ }} /usr/include /usr/include/hdf5/serial; do
      [ -d "$d" ] && INCDIRS="$INCDIRS -I$d"
    done
    HDF5LIBS="-lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5"
    FFTH=""
    for d in ${CPATH:+${CPATH//:/ }} /usr/include; do
      if [ -f "$d/fftw3.f03" ]; then FFTH="$d/fftw3.f03"; break; fi
    done
    MPIH=""
    ;;
esac
[ -n "$FC" ] && [ -x "$FC" ] || fail "no Fortran compiler found (toolchain $TOOLCHAIN)"

# --- MPI works? (two-task hello world) ----------------------------------------
MPIENV=""
if [ "$WITH_MPI" -eq 1 ] && want build; then
  if [ -z "$MPIFC" ] || [ ! -x "$MPIFC" ] || [ -z "$MPIRUN" ]; then
    warn "no MPI found: building the serial version only"; WITH_MPI=0
  else
    T="$(mktemp -d)"
    printf 'program t\nuse mpi\ninteger :: ierr, r\ncall MPI_Init(ierr)\ncall MPI_Comm_rank(MPI_COMM_WORLD,r,ierr)\nprint *,"VESPAWRANK",r\ncall MPI_Finalize(ierr)\nend program t\n' > "$T/t.f90"
    if ( cd "$T" && OMPI_FC="$FC" "$MPIFC" -o t.x t.f90 ) >>"$LOG" 2>&1 \
       && [ "$(cd "$T" && OMPI_MCA_rmaps_base_oversubscribe=1 perl -e 'alarm 60; exec @ARGV' "$MPIRUN" -np 2 ./t.x 2>>"$LOG" | grep -c VESPAWRANK)" = 2 ]; then
      say "MPI works ($MPIRUN)"
    else
      warn "MPI does not run on this machine ($MPIRUN: two-task test failed, see $LOG)."
      warn "Building the serial version only. Fix MPI and rerun ./install.sh to add it."
      WITH_MPI=0
    fi
    rm -rf "$T"
  fi
fi

# -----------------------------------------------------------------------------
#  3. parameter file and build
# -----------------------------------------------------------------------------
ARCH="gfortran_${PLAT//-/_}"
PARMS="$REPO/parms.auto"
TPAR="F"; [ "$WITH_MPI" -eq 1 ] && TPAR="T"
cat > "$PARMS" <<EOF
# generated by install.sh on $(date '+%Y-%m-%d %H:%M') for $PLAT ($TOOLCHAIN) -- rerun install.sh to change
ARCH="$ARCH"
TPARALLEL="$TPAR"
SPECIAL="none"
BLASDIR="$TOOL"
LAPACKDIR="$TOOL"
FFTDIR="$TOOL"
MPIDIR="$TOOL"
FFT_HEADER="$FFTH"
MPI_HEADER="$MPIH"
COMPILER_SCALAR="$FC -DCPPVAR_COMPILER_G95 $INCDIRS -fopenmp"
COMPILER_PARALLEL="OMPI_FC=$FC $MPIFC -DCPPVAR_COMPILER_G95 $INCDIRS -fopenmp"
FCFLAGS_NONE="-c "
FCFLAGS_OPT="-c -O3 -march=native -funroll-loops -fomit-frame-pointer -fimplicit-none -std=legacy"
FCFLAGS_PROF="-c -pg -O3 -funroll-loops"
FCFLAGS_DBG="-c -g -fbounds-check -fimplicit-none -std=legacy"
LDFLAGS_SCALAR="-I\${OBJDIR} $INCDIRS -L\${OBJDIR} $LIBDIRS"
LDFLAGS_PARALLEL="-I\${OBJDIR} $INCDIRS -L\${OBJDIR} $LIBDIRS"
LIBS_SCALAR="-lfftw3 -lopenblas -lpthread -lm -fopenmp $HDF5LIBS"
LIBS_PARALLEL="-lfftw3 -lopenblas -lpthread -lm -fopenmp $HDF5LIBS"
CPPFLAGS="-DCPPVAR_FFT_FFTW3 -DCPPVAR_BLAS_OPENBLAS -DCPPVAR_HDF5 -DCPPVAR_U77 -DCPPVAR_COMPILER_G95 "
FEXT="f90"
EOF
BINDIR="$REPO/bin/$ARCH"

if want build; then
  say "configuring (parameter file parms.auto)"
  export PATH="$TOOL/bin:$PATH"
  ( cd "$REPO" && ./configure --with-parmfile="$PARMS" ) >>"$LOG" 2>&1 || fail "configure failed"
  chmod +x "$BINDIR"/f90pp 2>/dev/null || true
  say "building the serial version (about 5-10 minutes)"
  ( cd "$REPO" && make fast_new ) >>"$LOG" 2>&1 || fail "serial build failed"
  if [ "$WITH_MPI" -eq 1 ]; then
    say "building the MPI version"
    ( cd "$REPO" && make fast_parallel_new ) >>"$LOG" 2>&1 || fail "MPI build failed"
  fi
  ( cd "$REPO" && make scripts ) >>"$LOG" 2>&1 || fail "installing the scripts failed"
fi

# -----------------------------------------------------------------------------
#  4. launcher
# -----------------------------------------------------------------------------
if want install; then
  cat > "$PREFIX/bin/vespaw" <<EOF
#!/usr/bin/env bash
# VESPAW launcher written by install.sh ($TOOLCHAIN toolchain): makes the
# build, its MPI and the setup library visible, then calls the dispatcher.
export VESPAW_HOME="$REPO"
export VESPAW_BIN_DIR="$BINDIR"
export VESPAW_POT_DIR="\${VESPAW_POT_DIR:-$REPO/vespaw_pot}"
export VESPAW_MPIRUN="\${VESPAW_MPIRUN:-$MPIRUN}"
export PATH="$BINDIR:$(dirname "${MPIRUN:-$FC}"):$ENV/bin:\$PATH"
export PYTHONPATH="$REPO/python\${PYTHONPATH:+:\$PYTHONPATH}"
if [ "\${1:-}" = "doctor" ]; then shift; exec "$REPO/installation/vespaw_doctor.sh" "\$@"; fi
exec "$BINDIR/vespaw" "\$@"
EOF
  chmod +x "$PREFIX/bin/vespaw"
  say "installed the launcher $PREFIX/bin/vespaw"
fi

# -----------------------------------------------------------------------------
#  5. check
# -----------------------------------------------------------------------------
if want test && [ "$RUN_TESTS" -eq 1 ]; then
  say "checking the build against the regression suite (a few minutes)"
  export VESPAW_BIN_DIR="$BINDIR"
  CHK="$PREFIX/check.log"
  ( cd "$REPO/regression" && "$PY" run.py pbe_atom si8_gamma ) >"$CHK" 2>&1 || true
  grep -E "\[FAIL\]|regression cases passed" "$CHK" || true
  grep -q "2/2 regression cases passed" "$CHK" || fail "serial check failed (see $CHK)"
  if [ "$WITH_MPI" -eq 1 ] && [ -x "$BINDIR/vespaw_fast_mpi.x" ]; then
    say "checking the MPI version on 2 tasks"
    ( cd "$REPO/regression" && VESPAW_REGRESSION_NP=2 VESPAW_MPIRUN="$MPIRUN" \
        "$PY" run.py si8_gamma ) >"$CHK.mpi" 2>&1 || true
    grep -E "\[FAIL\]|regression cases passed" "$CHK.mpi" || true
    grep -q "1/1 regression cases passed" "$CHK.mpi" || fail "MPI check failed (see $CHK.mpi)"
  fi
fi

say "done."
cat <<EOF

  Add VESPAW to your PATH (e.g. in ~/.zshrc or ~/.bashrc):
      export PATH="$PREFIX/bin:\$PATH"
  then:
      vespaw doctor          check the installation
      vespaw FILE.cntl       run a calculation
  Full log: $LOG
EOF
