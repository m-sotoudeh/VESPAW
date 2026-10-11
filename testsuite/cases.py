"""Multi-run regression cases driven by a case.json file.

A case directory contains template inputs (e.g. si.cntl, si.strc) and a
case.json:

{
  "description": "...",
  "inputs": ["si.cntl", "si.strc"],      # first entry is the .cntl to run
  "timeout": 900,                        # seconds per stage (optional)
  "runs": {
    "dav_var": {"subs": {"@PSIDYN@": "OPTIMIZER='DAV' MIXER='VARIATIONAL'"}},
    "nio_seeded": {"stages": [ {"subs": {...}}, {"subs": {...}} ]}
  },
  "defaults": {"@START@": "T", ...},     # substitutions applied to every run
  "checks": [
    {"run": "dav_var", "quantity": "energy", "value": -8.0070292, "tol": 2e-6,
     "reference": "where the value comes from"}
  ],
  "agree": [
    {"runs": ["dav_var", "rmm_var"], "quantity": "energy", "tol": 2e-6}
  ],
  "known_issues": {"dav_var_tscan": "why this run is expected to misbehave"},
  "fd": [
    {"quantity": "force:SI2:x", "run": "c", "plus": "xp", "minus": "xm",
     "step": 0.0188973, "factor": -1.0, "tol": 5e-6}
  ]
}

Quantities: "energy" (final E in Ha from the finish line), "gap_ev" (last
"ABSOLUTE GAP" in the .prot), "iterations" (from the finish line), and,
for runs with !CONTROL!GENERIC FORCES=T / STRESS=T, "force:<ATOM>:<x|y|z>"
[H/bohr], "dedeps:<i><j>" and "dedeps_trace" (dE/depsilon [H]) and
"pressure_gpa".

An "fd" entry compares an analytic quantity of run "run" with the central
finite difference factor*(E[plus]-E[minus])/(2*step) of two other runs.

A run with "optional": true (e.g. the long runs behind the tables of a
paper) is skipped unless it is named explicitly ("case:run") or the
environment variable VESPAW_REGRESSION_OPTIONAL=1 is set.

Runs listed in "known_issues" are executed and reported, but their checks
do not count towards the pass/fail result. All stages of one run share a
work directory, so later stages can continue from the restart file of an
earlier one (use "START=F" in their substitutions).
"""
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

# VESPAW prints energies in eV, forces in eV/A and dE/depsilon in eV since
# 2026-10-02; the references here stay in atomic units (CODATA 2018 factors).
_HARTREE_EV = 27.211386245988
_BOHR_A = 0.529177210903
_FINISH_RE = re.compile(r"AFTER\s+(\d+)\s+ITERATIONS\s*\|\s*E\s*=\s*([+-]?\d+\.\d+)\s*EV")
_GAP_RE = re.compile(r"ABSOLUTE GAP\.*:\s*([+-]?\d+\.\d+)\s*EV")
_FORCE_RE = re.compile(r"^FORCE (\S+)\s+(\S+)\s+(\S+)\s+(\S+)", re.M)
_STRESS_RE = re.compile(r"^STRESS \s*\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)", re.M)
_PRESSURE_RE = re.compile(r"^PRESSURE\[GPA\]\s+(\S+)", re.M)


def _substitute(text, subs):
    for key, val in subs.items():
        text = text.replace(key, val)
    return text


def _run_stage(case_dir, work, inputs, subs, env, timeout, binary):
    for name in inputs:
        text = (case_dir / name).read_text()
        (work / name).write_text(_substitute(text, subs))
    unresolved = [n for n in inputs if re.search(r"@[A-Z_]+@", (work / n).read_text())]
    if unresolved:
        raise RuntimeError(f"unresolved @PLACEHOLDER@ in {', '.join(unresolved)}")
    cntl = inputs[0]
    t0 = time.time()
    cmd = [binary, cntl]
    nproc = int(os.environ.get("VESPAW_REGRESSION_NP", "1"))
    if nproc > 1:
        # MPI run: the MPI executable next to the serial one, one thread per rank
        mpibin = str(Path(binary).with_name(Path(binary).stem + "_mpi.x"))
        mpirun = (os.environ.get("VESPAW_MPIRUN") or shutil.which("mpirun")
                  or "/opt/homebrew/bin/mpirun")
        cmd = [mpirun, "--oversubscribe", "-np", str(nproc), mpibin, cntl]
        env = dict(env, OMP_NUM_THREADS="1")
    proc = subprocess.run(cmd, cwd=work, env=env, capture_output=True,
                          text=True, timeout=timeout)
    wall = time.time() - t0
    (work / "run.log").write_text(proc.stdout + proc.stderr)
    return proc, wall


def _parse(work, stdout):
    out = {}
    matches = _FINISH_RE.findall(stdout)
    if matches:
        out["iterations"] = int(matches[-1][0])
        out["energy"] = float(matches[-1][1]) / _HARTREE_EV
    for name, fx, fy, fz in _FORCE_RE.findall(stdout):
        for comp, val in zip("xyz", (fx, fy, fz)):
            out[f"force:{name}:{comp}"] = float(val) * _BOHR_A / _HARTREE_EV
    rows = _STRESS_RE.findall(stdout)
    if len(rows) >= 3:
        rows = rows[-3:]
        for i in range(3):
            for j in range(3):
                out[f"dedeps:{i + 1}{j + 1}"] = float(rows[i][j]) / _HARTREE_EV
        out["dedeps_trace"] = sum(float(rows[i][i]) for i in range(3)) / _HARTREE_EV
    pres = _PRESSURE_RE.findall(stdout)
    if pres:
        out["pressure_gpa"] = float(pres[-1])
    # the detailed legacy protocol (with the ABSOLUTE GAP line) is the .log file
    prot = next((q for q in work.glob("*.log") if q.name != "run.log"), None)
    if prot is not None:
        gaps = _GAP_RE.findall(prot.read_text(errors="replace"))
        if gaps:
            out["gap_ev"] = float(gaps[-1])
    return out


def run_case(case_dir, env, binary, only_runs=None, keep=False):
    spec = json.loads((case_dir / "case.json").read_text())
    inputs = spec["inputs"]
    timeout = spec.get("timeout", 900)
    defaults = spec.get("defaults", {})
    known = spec.get("known_issues", {})
    results = {}
    print(f"--- {case_dir.name}: {spec.get('description', '')}")

    for run_name, run in spec["runs"].items():
        if only_runs and run_name not in only_runs:
            continue
        if (run.get("optional") and not only_runs
                and os.environ.get("VESPAW_REGRESSION_OPTIONAL") != "1"):
            print(f"  [SKIP] {run_name}: optional run (VESPAW_REGRESSION_OPTIONAL=1 to include)")
            continue
        stages = run.get("stages", [run])
        work = case_dir / "work" / run_name
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        total = 0.0
        res = {}
        failed = None
        for i, stage in enumerate(stages):
            subs = dict(defaults)
            subs.update(stage.get("subs", {}))
            try:
                proc, wall = _run_stage(case_dir, work, inputs, subs, env, timeout, binary)
            except subprocess.TimeoutExpired:
                failed = f"stage {i + 1} timed out after {timeout} s"
                break
            total += wall
            if proc.returncode != 0:
                failed = f"stage {i + 1} exited {proc.returncode}"
                break
            res = _parse(work, proc.stdout)
            if "energy" not in res:
                failed = f"stage {i + 1}: no finish line"
                break
        res["wall_s"] = total
        res["failed"] = failed
        results[run_name] = res
        tag = " (known issue)" if run_name in known else ""
        if failed:
            print(f"    {run_name}: FAILED TO RUN: {failed}{tag}")
        else:
            gap = f", gap {res['gap_ev']:.4f} eV" if "gap_ev" in res else ""
            print(f"    {run_name}: E = {res['energy'] * _HARTREE_EV:.8f} eV, {res.get('iterations', '?')} it"
                  f"{gap}, {total:.0f} s{tag}")
        if not keep:
            for f in work.glob("*.rstrt"):
                f.unlink()
            for f in work.glob("*.tra"):
                f.unlink()

    ok_all = True
    for chk in spec.get("checks", []):
        run_name = chk["run"]
        if run_name not in results:
            continue
        res = results[run_name]
        q = chk["quantity"]
        counts = run_name not in known
        if res["failed"] or q not in res:
            status = "FAIL"
            detail = f"{q} not available ({res['failed'] or 'not in output'})"
            ok = False
        else:
            diff = abs(res[q] - chk["value"])
            ok = diff <= chk["tol"]
            status = "PASS" if ok else "FAIL"
            detail = (f"{q} = {res[q]:.10g} (expected {chk['value']:.10g} +/- "
                      f"{chk['tol']:g}, diff {diff:.3g})")
        if not counts:
            status = "KNOWN-" + status
        print(f"  [{status}] {run_name}: {detail}")
        if counts and not ok:
            ok_all = False

    for agr in spec.get("agree", []):
        runs = [r for r in agr["runs"] if r in results]
        if len(runs) < 2:
            continue
        q = agr["quantity"]
        vals = {r: results[r].get(q) for r in runs if not results[r]["failed"]}
        counts = not any(r in known for r in runs)
        if len(vals) < len(runs) or any(v is None for v in vals.values()):
            ok = False
            detail = f"{q} missing for some of {', '.join(runs)}"
        else:
            spread = max(vals.values()) - min(vals.values())
            ok = spread <= agr["tol"]
            detail = f"{q} spread over {', '.join(runs)} = {spread:.3g} (tol {agr['tol']:g})"
        status = "PASS" if ok else "FAIL"
        if not counts:
            status = "KNOWN-" + status
        print(f"  [{status}] agreement: {detail}")
        if counts and not ok:
            ok_all = False

    for fd in spec.get("fd", []):
        need = [fd["run"], fd["plus"], fd["minus"]]
        if not all(r in results for r in need):
            continue
        q = fd["quantity"]
        bad = [r for r in need if results[r]["failed"]]
        if bad or q not in results[fd["run"]]:
            ok = False
            detail = f"{q}: missing ({', '.join(bad) or 'not in output'})"
        else:
            num = fd.get("factor", 1.0) * (results[fd["plus"]]["energy"]
                                           - results[fd["minus"]]["energy"]) / (2.0 * fd["step"])
            ana = results[fd["run"]][q]
            ok = abs(ana - num) <= fd["tol"]
            detail = (f"{q} = {ana:.8g} vs finite difference {num:.8g} "
                      f"(diff {abs(ana - num):.3g}, tol {fd['tol']:g})")
        print(f"  [{'PASS' if ok else 'FAIL'}] fd: {detail}")
        if not ok:
            ok_all = False

    for run_name, why in known.items():
        if run_name in results:
            print(f"  note: {run_name} is a known issue: {why}")
    return ok_all
