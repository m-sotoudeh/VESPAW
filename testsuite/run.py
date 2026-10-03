#!/usr/bin/env python3
"""VESPAW regression guard -- a small, fast, checked-in set of known-good
energies, distinct from the (gitignored) exploratory `tests/` tree.

Why this exists: every session's "verified against a real run" was ad hoc
and left no trace in the repository (tests/ is gitignored). This is the
minimum needed so a future Phase 0.5 cleanup pass, or any other change,
has something automatic to check against instead of relying on someone
remembering to re-run last session's manual checks. See the dev notes's
"Not yet implemented -- continuation checklist" for the YAML frontend, and
the "Known bugs" section for SCAN0r, for why these particular cases were
picked.

Each case directory holds its own real input (.cntl/.strc, or a YAML
input) plus `expected.json` (a known-good final energy + tolerance, in
Hartree). Run-generated files (.rstrt/.prot/results.h5/etc.) are not
committed -- see testsuite/.gitignore.

Usage:
    python3 testsuite/run.py                 # run every case
    python3 testsuite/run.py pbe_atom         # run just one case
    python3 testsuite/run.py si_pbe:dav_var   # one run of a multi-run case
    python3 testsuite/run.py --list           # list all cases

Multi-run cases (a directory with case.json, see cases.py) solve the same
system with several optimizers and check energies, gaps and cross-method
agreement against references with recorded provenance.

Requires the `vespaw` dispatcher on PATH and VESPAW_POT_DIR pointing at
../potentials (both set automatically here, relative to this file).
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

REGRESSION_DIR = Path(__file__).resolve().parent
REPO_ROOT = REGRESSION_DIR.parent


def _find_bin_dir():
    """The build directory bin/<arch>/ holding vespaw_fast.x: $VESPAW_BIN_DIR if
    set, otherwise the most recently built one."""
    env = os.environ.get("VESPAW_BIN_DIR")
    if env:
        return Path(env)
    found = sorted((REPO_ROOT / "bin").glob("*/vespaw_fast.x"), key=lambda p: p.stat().st_mtime)
    return found[-1].parent if found else REPO_ROOT / "bin" / "gfortran_openblas_fftw3_arm64"


VESPAW_BIN_DIR = _find_bin_dir()
VESPAW_POT_DIR = REPO_ROOT / "potentials"
PYTHON_PKG_DIR = REPO_ROOT / "python"

# VESPAW prints and stores energies in eV since 2026-10-02; the expected
# values here stay in Hartree (CODATA 2018 factor).
_HARTREE_EV = 27.211386245988
_FINISH_RE = re.compile(r"ITERATIONS\s*\|\s*E\s*=\s*([+-]?\d+\.\d+)\s*EV\b")


def _run_env():
    env = dict(os.environ)
    env["PATH"] = f"{VESPAW_BIN_DIR}{os.pathsep}{env.get('PATH', '')}"
    env["VESPAW_POT_DIR"] = str(VESPAW_POT_DIR)
    pypath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{PYTHON_PKG_DIR}{os.pathsep}{pypath}" if pypath else str(PYTHON_PKG_DIR)
    return env


def _check(name, energy, expected):
    tol = expected["tolerance_ha"]
    want = expected["energy_ha"]
    diff = abs(energy - want)
    ok = diff <= tol
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}: E = {energy:.10f} Ha "
          f"(expected {want:.10f} +/- {tol:g}, diff {diff:.3g})")
    return ok


def run_cntl_case(case_dir):
    """Cases driven directly by a .cntl deck (found by glob, one per dir)."""
    expected = json.loads((case_dir / "expected.json").read_text())
    cntl = next(case_dir.glob("*.cntl"))
    proc = subprocess.run(["vespaw", cntl.name], cwd=case_dir, env=_run_env(),
                           capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        print(f"[FAIL] {case_dir.name}: vespaw exited {proc.returncode}")
        print(proc.stdout[-2000:])
        print(proc.stderr[-2000:])
        return False
    m = _FINISH_RE.search(proc.stdout)
    if not m:
        print(f"[FAIL] {case_dir.name}: no 'E = ... EV' finish line found in output")
        return False
    return _check(case_dir.name, float(m.group(1)) / _HARTREE_EV, expected)


def run_yaml_frontend_case(case_dir):
    sys.path.insert(0, str(PYTHON_PKG_DIR))
    from vespaw.yaml_input import prepare_run  # noqa: E402  (path set above)

    expected = json.loads((case_dir / "expected.json").read_text())
    yaml_path = next(case_dir.glob("*.yaml"))
    env = _run_env()
    old_environ = dict(os.environ)
    os.environ.update(env)
    try:
        calc = prepare_run(yaml_path)
        result = calc.run(timeout=300, write_run_json=False)
    finally:
        os.environ.clear()
        os.environ.update(old_environ)
    return _check(case_dir.name, result.energy / _HARTREE_EV, expected)


CASES = {
    "pbe_atom": run_cntl_case,
    "pbesl0r_atom": run_cntl_case,
    "yaml_frontend": run_yaml_frontend_case,
    "yaml_raw": run_yaml_frontend_case,
}


def _multi_run_cases():
    """Case directories driven by case.json (see cases.py)."""
    return sorted(p.parent.name for p in REGRESSION_DIR.glob("*/case.json"))


def run_multi_case(case_dir, only_runs=None):
    from cases import run_case  # noqa: E402  (local module)
    env = _run_env()
    env.setdefault("OMP_NUM_THREADS", "8")
    return run_case(case_dir, env, str(VESPAW_BIN_DIR / "vespaw_fast.x"),
                    only_runs=only_runs)


def main(argv):
    for name in _multi_run_cases():
        CASES.setdefault(name, run_multi_case)
    if argv and argv[0] == "--list":
        for name in CASES:
            print(name)
        return 0
    # NAME:RUN selects a single run of a multi-run case
    selected = {}
    for arg in argv:
        name, _, run = arg.partition(":")
        selected.setdefault(name, set())
        if run:
            selected[name].add(run)
    names = list(selected) or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        print(f"unknown case(s): {', '.join(unknown)} (known: {', '.join(CASES)})")
        return 1

    results = {}
    for name in names:
        case_dir = REGRESSION_DIR / name
        try:
            if CASES[name] is run_multi_case:
                results[name] = run_multi_case(case_dir, selected.get(name) or None)
            else:
                results[name] = CASES[name](case_dir)
        except Exception as exc:  # noqa: BLE001 -- report, don't crash the suite
            print(f"[FAIL] {name}: raised {exc!r}")
            results[name] = False

    print()
    n_pass = sum(results.values())
    print(f"{n_pass}/{len(results)} regression cases passed")
    return 0 if n_pass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
