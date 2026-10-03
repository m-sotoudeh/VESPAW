#!/usr/bin/env bash
# =============================================================================
#  vespaw doctor -- report the state of a VESPAW installation.
#  Called by the launcher written by install.sh ("vespaw doctor"); prints
#  what was found and a one-line verdict. --run also runs a 1-minute check.
# =============================================================================
REPO="${VESPAW_HOME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BIN="${VESPAW_BIN_DIR:-}"
if [ -z "$BIN" ]; then
  BIN="$(ls -td "$REPO"/bin/*/ 2>/dev/null | head -1)"; BIN="${BIN%/}"
fi
ok=0; warn=0; bad=0
pass() { printf '  \033[32mok\033[0m    %s\n' "$*"; ok=$((ok+1)); }
note() { printf '  \033[33mwarn\033[0m  %s\n' "$*"; warn=$((warn+1)); }
miss() { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; bad=$((bad+1)); }

echo "VESPAW doctor"
echo "  source      $REPO"
echo "  build       ${BIN:-none}"
echo "  platform    $(uname -s) $(uname -m)"
echo

# --- executables -------------------------------------------------------------
if [ -x "$BIN/vespaw_fast.x" ]; then pass "serial executable vespaw_fast.x"
else miss "serial executable missing (run ./install.sh)"; fi
if [ -x "$BIN/vespaw_fast_mpi.x" ]; then pass "MPI executable vespaw_fast_mpi.x"
else note "no MPI executable (serial runs only)"; fi
if [ -x "$BIN/vespaw" ]; then pass "dispatcher bin/.../vespaw"
else miss "dispatcher missing (make scripts)"; fi

# --- shared libraries ----------------------------------------------------------
if [ -x "$BIN/vespaw_fast.x" ]; then
  if [ "$(uname -s)" = "Darwin" ]; then libs="$(otool -L "$BIN/vespaw_fast.x" 2>/dev/null)"
  else libs="$(ldd "$BIN/vespaw_fast.x" 2>/dev/null)"; fi
  if echo "$libs" | grep -q "not found"; then
    miss "unresolved libraries:"; echo "$libs" | grep "not found" | sed 's/^/          /'
  else pass "all shared libraries resolve"; fi
  for l in openblas fftw3 hdf5 gomp gfortran; do
    echo "$libs" | grep -q "$l" || note "library $l not linked"
  done
  if [ "$(uname -s)" = "Darwin" ]; then
    a="$(lipo -archs "$BIN/vespaw_fast.x" 2>/dev/null)"
    if [ "$(uname -m)" = "arm64" ] && [ "$a" != "arm64" ]; then
      note "executable is $a on an arm64 Mac (runs under Rosetta, ~2x slower)"
    else pass "native architecture ($a)"; fi
  fi
fi

# --- MPI -------------------------------------------------------------------------
if command -v mpirun >/dev/null 2>&1; then
  MR="$(command -v mpirun)"
  if [ "$(uname -s)" = "Darwin" ]; then
    MRREAL="$(python3 -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' "$MR" 2>/dev/null || echo "$MR")"
    MA="$(lipo -archs "$MRREAL" 2>/dev/null)"
    if [ -n "$MA" ] && [ "$MA" != "$(uname -m)" ]; then
      note "mpirun ($MR) is built for $MA, not $(uname -m): it will not launch the native build; put the matching MPI first in PATH"
    else pass "mpirun: $MR"; fi
  else pass "mpirun: $MR"; fi
else note "mpirun not on PATH"; fi

# --- setups, python ----------------------------------------------------------------
POT="${VESPAW_POT_DIR:-$REPO/potentials}"
n=$(ls "$POT"/*.strc-species 2>/dev/null | wc -l | tr -d ' ')
if [ "$n" -gt 0 ]; then pass "setup library: $n species files in $POT"
else miss "setup library not found ($POT)"; fi
if python3 -c "import numpy, yaml" 2>/dev/null; then pass "python3 with numpy and pyyaml"
else note "python3 lacks numpy/pyyaml (YAML input and Python API unavailable)"; fi

# --- threads -------------------------------------------------------------------------
echo
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}  OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-unset}"
echo "  (with MPI, set both to 1 unless you deliberately combine MPI and threads)"

# --- optional quick run ----------------------------------------------------------------
if [ "${1:-}" = "--run" ] && [ -x "$BIN/vespaw_fast.x" ]; then
  echo; echo "  running a 1-minute check (regression case pbe_atom)..."
  if ( cd "$REPO/testsuite" && VESPAW_BIN_DIR="$BIN" python3 run.py pbe_atom ) 2>&1 \
       | grep -q "1/1 regression cases passed"; then pass "check calculation"
  else miss "check calculation failed (python3 testsuite/run.py pbe_atom)"; fi
fi

echo
if [ "$bad" -gt 0 ]; then echo "  verdict: NOT working ($bad problem(s)); see above"; exit 1
elif [ "$warn" -gt 0 ]; then echo "  verdict: working, with $warn warning(s)"
else echo "  verdict: working"; fi
