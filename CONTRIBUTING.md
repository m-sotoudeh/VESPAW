# Contributing to VESPAW

Bug reports, questions, test cases and code are welcome.

## Reporting a problem

Open an [issue](https://github.com/m-sotoudeh/VESPAW/issues) with the
input files (`.cntl`, `.strc` or YAML), the relevant part of the `.prot`
report and the output of `vespaw doctor`. A small system that shows the
problem is the most useful. Known problems are listed in
[doc/KNOWN_ISSUES.md](doc/KNOWN_ISSUES.md).

## Changing the code

1. Build with `./install.sh` (see [doc/INSTALL.md](doc/INSTALL.md)).
2. Make the change on a branch. Follow the style of the surrounding code:
   VESPAW keeps the CP-PAW conventions (upper-case free-form Fortran,
   object-style modules with `$`-separated interfaces such as
   `WAVES$SETR8`). OpenMP directives must be written in lower case
   (`!$omp ...`, including continuation lines); the preprocessor keeps
   only those.
3. Run the regression suite:

   ```sh
   python3 testsuite/run.py            # all cases
   python3 testsuite/run.py si_pbe     # one case
   ```

   Every case solves the same system with several optimizers and checks
   energies, gaps, forces and stresses against references. If a reference
   changes on purpose, update `case.json` and state in its `reference`
   field what changed and why.
4. Document new input keywords in `doc/src/manual` and add a line to
   `CHANGELOG.md` under "Unreleased".
5. Open a pull request describing the change and how it was checked.

## New functionals and methods

Please include a test case: a small system with a reference value from an
independent source (another code, libxc, a finite-difference check or a
publication), so the implementation stays verified.

## License

Contributions are distributed under the GNU General Public License,
version 3, like the rest of VESPAW (see [LICENSE](LICENSE) and
[NOTICE](NOTICE)).
