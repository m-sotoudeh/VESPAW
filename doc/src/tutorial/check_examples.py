#!/usr/bin/env python3
"""Run the tutorial examples and check their results.

Each directory examples/NN_name/ is a case in the format of the
regression suite (see cases.py there): template inputs, runs and
expected values. Usage:

    python3 doc/tutorial/check_examples.py            # all examples
    python3 doc/tutorial/check_examples.py 02_si      # one example
    python3 doc/tutorial/check_examples.py 02_si:scf  # one run
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "regression").exists():          # public layout: doc/src/tutorial
    ROOT = Path(__file__).resolve().parents[3]
SUITE = ROOT / "regression" if (ROOT / "regression").exists() else ROOT / "testsuite"
sys.path.insert(0, str(SUITE))
from cases import run_case  # noqa: E402


def main(argv):
    examples = Path(__file__).resolve().parent / "examples"
    if not examples.exists():                     # public layout: examples/tutorial
        examples = ROOT / "examples" / "tutorial"
    built = sorted((ROOT / "bin").glob("*/vespaw_fast.x"), key=lambda p: p.stat().st_mtime)
    binary = os.environ.get("VESPAW_BINARY", str(built[-1]) if built else "vespaw_fast.x")
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "8")
    selected = {}
    for arg in argv:
        name, _, run = arg.partition(":")
        selected.setdefault(name, set())
        if run:
            selected[name].add(run)
    ok = True
    for case in sorted(p.parent for p in examples.glob("*/case.json")):
        if selected and case.name not in selected:
            continue
        ok &= run_case(case, env, binary, only_runs=selected.get(case.name) or None)
    print("ALL EXAMPLES OK" if ok else "SOME EXAMPLES FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
