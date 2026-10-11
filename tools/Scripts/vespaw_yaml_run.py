#!/usr/bin/env python3
"""vespaw run -- validate, lower, and run a YAML input file end to end.

Thin CLI over vespaw.yaml_input.prepare_run() + Calculation.run(): loads
input.yaml, validates it, lowers it to a .cntl/.strc pair next to the YAML
file (or under --workdir), runs the real vespaw binary with its output
streamed live to this terminal (so the live-SCF table shows as the run
progresses, not dumped all at once at the end), and prints a short summary
on success.

See python/vespaw/yaml_input.py's module docstring for the input schema
("simple case" fields plus the raw: passthrough that reaches every .cntl
keyword).
"""
import argparse
import sys
from pathlib import Path

# This script lives at <repo>/tools/Scripts/vespaw_yaml_run.py (and gets
# copied verbatim into <repo>/bin/<arch>/ by `make scripts`, so the layout
# it runs from is not fixed) -- climb from this file to the actual repo
# checkout rather than assuming a fixed relative depth.
def _find_repo_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "python" / "vespaw" / "__init__.py").exists():
            return candidate
        if (candidate / "src" / "Tools" / "Scripts" / "vespaw_convert.py").exists() \
                and (candidate.parent / "python" / "vespaw" / "__init__.py").exists():
            return candidate.parent
    raise SystemExit(
        "vespaw run: could not locate the VESPAW checkout's python/ "
        "directory from this script's location -- expected to find "
        "python/vespaw/__init__.py somewhere above "
        f"{start}")


_REPO_ROOT = _find_repo_root(Path(__file__).resolve().parent)
sys.path.insert(0, str(_REPO_ROOT / "python"))

from vespaw import VespawError  # noqa: E402
from vespaw.yaml_input import VespawInputError, prepare_run  # noqa: E402

_HARTREE_IN_EV = 27.211386245988


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="vespaw run",
        description="Validate, lower, and run a VESPAW YAML input file.")
    p.add_argument("input", help="path to input.yaml")
    p.add_argument("--workdir", default=None,
                    help="directory for the generated .cntl/.strc and "
                         "results (default: alongside the input file)")
    p.add_argument("--nprocs", type=int, default=1,
                    help="run under 'mpirun -np NPROCS' (default: 1, serial)")
    p.add_argument("--timeout", type=float, default=None,
                    help="wall-clock timeout in seconds (default: none)")
    p.add_argument("--no-run-json", action="store_true",
                    help="skip writing the run.json provenance record next "
                         "to the results")
    args = p.parse_args(argv)

    yaml_path = Path(args.input)
    if not yaml_path.exists():
        print(f"vespaw run: no such file: {yaml_path}", file=sys.stderr)
        return 1

    try:
        calc = prepare_run(
            yaml_path,
            workdir=Path(args.workdir) if args.workdir else None)
    except VespawInputError as exc:
        print(f"VESPAW INPUT ERROR\n\n{exc}\n\nCalculation not started.",
              file=sys.stderr)
        return 1

    # Calculation.run() otherwise looks up "vespaw" on $PATH (shutil.which)
    # and fails if the caller invoked this very script by its full path
    # without ever putting bin/<arch>/ on $PATH -- but the dispatcher is
    # sitting right next to this script (both land in bin/<arch>/ via
    # `make scripts`), so there's no real ambiguity to resolve here.
    colocated_vespaw = Path(__file__).resolve().parent / "vespaw"
    if colocated_vespaw.exists():
        calc.vespaw_bin = str(colocated_vespaw)

    calc.nprocs = args.nprocs

    try:
        result = calc.run(timeout=args.timeout, stream=True,
                           write_run_json=not args.no_run_json)
    except VespawError as exc:
        print(f"\nvespaw run: {exc}", file=sys.stderr)
        return 1

    print()
    print(f"Energy : {result.energy:.8f} eV")
    if result.gap is not None:
        print(f"Gap    : {result.gap:.4f} eV")
    elif result.is_metal:
        print("Gap    : metallic (no gap)")
    print(f"Results: {calc.hdf5}")
    print(f"Report : {calc.cntl.with_suffix('.prot')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
