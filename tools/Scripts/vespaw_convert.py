#!/usr/bin/env python3
"""vespaw convert -- build a VESPAW/CP-PAW .strc file's structural part
(!GENERIC/!OCCUPATIONS/!KPOINTS/!LATTICE/!ATOM) from any structure file
ASE can read (.cif, POSCAR/CONTCAR, .xyz, extxyz, ...), so positions and
lattice vectors never have to be hand-typed.

What this does NOT do, on purpose: invent !SPECIES (PAW setup) blocks.
Those encode real, validated pseudopotential physics (NDLSS/NTBO onsite-
hybrid parameters) that a structure file has no way to supply -- they are
copied verbatim, per chemical element, out of one or more existing .strc
files (or a potentials library file, see ../../../potentials/README.md) passed
via --species-from. If an element in the input structure has no matching
species block available, this tool refuses to guess and tells you exactly
which element is missing.

Element -> species-block matching is data-driven, not name-guessed: each
!SPECIES block is matched to a chemical element via its own !AUGMENT
sub-block, reading EL='<symbol>' when present, or Z=<atomic number>
otherwise (both conventions are used in the wild -- see e.g. the MG vs MN
!SPECIES blocks in tests/mgo/mgo.strc).
"""
import argparse
import re
import sys

TOKEN_RE = re.compile(r"'[^']*'|\S+")


def strip_comments(text):
    """Blank out '#'-led comment lines, exactly as the real .cntl/.strc
    parser treats them (src/vespaw_linkedlist.f90, BUFFER_READ) -- a
    comment must be the first non-blank character on its own line."""
    return '\n'.join('' if line.strip().startswith('#') else line
                      for line in text.splitlines())


def find_species_blocks(text):
    """Return the exact source text of every top-level '!SPECIES ... !END'
    block in `text`, found by depth-counting '!'-prefixed tokens ('!END'
    closes, anything else '!X' opens) -- the same bracketing rule the
    LINKEDLIST reader itself relies on, so this works regardless of
    whether the '!SPECIES' block sits inside a full !STRUCTURE or in a
    bare species-library file."""
    cleaned = strip_comments(text)
    toks = [(m.group(0), m.start(), m.end()) for m in TOKEN_RE.finditer(cleaned)]
    blocks = []
    i, n = 0, len(toks)
    while i < n:
        tok, start_off, end_off = toks[i]
        if tok.upper() == '!SPECIES':
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                t2, _s2, e2 = toks[j]
                if t2.upper() == '!END':
                    depth -= 1
                elif t2.startswith('!'):
                    depth += 1
                end_off = e2
                j += 1
            if depth != 0:
                raise ValueError("unterminated !SPECIES block (missing !END)")
            blocks.append(cleaned[start_off:end_off])
            i = j
        else:
            i += 1
    return blocks


def species_tag(block_text):
    m = re.search(r"\bNAME\s*=\s*'([^']*)'", block_text, re.IGNORECASE)
    if not m:
        raise ValueError("!SPECIES block has no NAME='...'")
    return m.group(1)


def species_element(block_text):
    """Chemical element for a !SPECIES block, from its !AUGMENT
    sub-block: EL='<symbol>' if present, else Z=<atomic number>."""
    m = re.search(r"\bEL\s*=\s*'([^']*)'", block_text, re.IGNORECASE)
    if m:
        return m.group(1).strip().upper()
    m = re.search(r"\bZ\s*=\s*([0-9.]+)", block_text, re.IGNORECASE)
    if m:
        from ase.data import chemical_symbols
        z = int(round(float(m.group(1))))
        if 0 < z < len(chemical_symbols):
            return chemical_symbols[z].upper()
    return None


def load_species_library(paths):
    """tag -> block_text and element -> tag, scanned across all
    --species-from files (later files may add elements the earlier ones
    lack; a duplicate element/tag combination is fine and just keeps the
    first one found)."""
    by_tag = {}
    by_element = {}
    for path in paths:
        with open(path) as fh:
            text = fh.read()
        for block in find_species_blocks(text):
            tag = species_tag(block)
            element = species_element(block)
            by_tag.setdefault(tag, block)
            if element and element not in by_element:
                by_element[element] = tag
    return by_tag, by_element


def build_strc(atoms, by_tag, by_element, element_map, args):
    symbols = atoms.get_chemical_symbols()
    unique_elements = sorted(set(s.upper() for s in symbols))

    missing = [el for el in unique_elements
               if el not in element_map and el not in by_element]
    if missing:
        available = ', '.join(sorted(by_element)) or '(none)'
        raise SystemExit(
            "vespaw convert: no !SPECIES block found for element(s) "
            f"{', '.join(missing)}.\n"
            f"  --species-from provided blocks for: {available}\n"
            "  Add a validated !SPECIES block for the missing element(s) "
            "(see potentials/README.md) and pass its file via --species-from, "
            "or map it explicitly with --map SYMBOL=TAG.")

    tag_of_element = dict(by_element)
    tag_of_element.update(element_map)  # --map always wins

    used_tags = []
    for el in unique_elements:
        tag = tag_of_element[el]
        if tag not in by_tag:
            raise SystemExit(
                f"vespaw convert: --map assigned element {el} to species "
                f"tag '{tag}', but no !SPECIES NAME='{tag}' block was found "
                "in any --species-from file.")
        used_tags.append((el, tag))

    species_text = '\n'.join('  ' + by_tag[tag] for _el, tag in used_tags)

    cell = atoms.get_cell()
    if not cell.any():
        raise SystemExit(
            "vespaw convert: input structure has no lattice vectors "
            "(cell is all-zero) -- an isolated/molecular structure needs "
            "an explicit vacuum cell before conversion (see ASE's "
            "`atoms.center(vacuum=...)`).")
    t_rows = [' '.join(f'{c:.10f}' for c in row) for row in cell]
    lattice_text = ('!LATTICE T= ' + t_rows[0] + '\n'
                     '              ' + t_rows[1] + '\n'
                     '              ' + t_rows[2] + ' !END')

    if args.kdiv is not None:
        kpoint_kw = 'DIV={} {} {}'.format(*args.kdiv)
    else:
        kpoint_kw = f'R={args.kr}'
    kpoint_kw += ' SHIFT={} {} {}'.format(*args.kshift)

    positions = atoms.get_positions()
    counts = {}
    atom_lines = []
    for sym, pos in zip(symbols, positions):
        tag = tag_of_element[sym.upper()]
        counts[tag] = counts.get(tag, 0) + 1
        name = f'{tag}{sum(counts.values())}'
        atom_lines.append(
            "  !ATOM NAME='{}' R= {:.10f} {:.10f} {:.10f} !END".format(name, *pos))

    lines = [
        '!STRUCTURE',
        '  !GENERIC LUNIT[AA]=1.0 !END',
        f'  !OCCUPATIONS EMPTY={args.empty} NSPIN={args.nspin} '
        f'SPIN[HBAR]={args.spin} CHARGE[E]={args.charge} !END',
        f'  !KPOINTS {kpoint_kw} !END',
        '  ' + lattice_text.replace('\n', '\n  '),
        species_text,
        '\n'.join(atom_lines),
        '!END',
        '!EOB',
    ]
    return '\n'.join(lines) + '\n', unique_elements, used_tags


def parse_map(pairs):
    out = {}
    for item in pairs or []:
        if '=' not in item:
            raise SystemExit(f"vespaw convert: --map expects SYMBOL=TAG, got '{item}'")
        sym, tag = item.split('=', 1)
        out[sym.strip().upper()] = tag.strip()
    return out


def main(argv=None):
    p = argparse.ArgumentParser(
        prog='vespaw convert',
        description="Build a VESPAW .strc file's structural blocks from a "
                    "CIF/POSCAR/xyz/... file, splicing in existing "
                    "!SPECIES (PAW setup) blocks by chemical element.")
    p.add_argument('input', help='structure file (any ASE-readable format)')
    p.add_argument('output', help='.strc file to write')
    p.add_argument('--species-from', nargs='+', required=True, metavar='FILE',
                    help='one or more .strc/species-library files to pull '
                         'validated !SPECIES blocks from')
    p.add_argument('--format', default=None,
                    help='force the ASE input format (default: guess from extension)')
    p.add_argument('--map', nargs='+', metavar='SYMBOL=TAG',
                    help="force element SYMBOL to species tag TAG (overrides "
                         "the automatic EL=/Z= match); use when a species "
                         "block has neither EL= nor a recognizable Z=")
    p.add_argument('--empty', type=int, default=4,
                    help='!OCCUPATIONS EMPTY= extra empty bands (default: 4)')
    p.add_argument('--nspin', type=int, default=1, choices=(1, 2, 3),
                    help='!OCCUPATIONS NSPIN= (default: 1)')
    p.add_argument('--spin', type=float, default=0.0,
                    help='!OCCUPATIONS SPIN[HBAR]= (default: 0.0)')
    p.add_argument('--charge', type=float, default=0.0,
                    help='!OCCUPATIONS CHARGE[E]= net charge in electrons (default: 0.0)')
    kpt = p.add_mutually_exclusive_group()
    kpt.add_argument('--kdiv', type=int, nargs=3, metavar=('NA', 'NB', 'NC'),
                      help='!KPOINTS DIV=na nb nc (default if neither --kdiv nor --kr given: 1 1 1)')
    kpt.add_argument('--kr', type=float, metavar='R',
                      help='!KPOINTS R= reciprocal-space radius (mutually exclusive with --kdiv)')
    p.add_argument('--kshift', type=int, nargs=3, default=(0, 0, 0),
                    metavar=('SA', 'SB', 'SC'), help='!KPOINTS SHIFT= (default: 0 0 0)')
    args = p.parse_args(argv)
    if args.kdiv is None and args.kr is None:
        args.kdiv = (1, 1, 1)

    try:
        from ase.io import read
    except ImportError:
        raise SystemExit(
            "vespaw convert: needs ASE (pip install ase) to read structure files")

    atoms = read(args.input, format=args.format)
    by_tag, by_element = load_species_library(args.species_from)
    element_map = parse_map(args.map)

    text, elements, used_tags = build_strc(atoms, by_tag, by_element, element_map, args)
    with open(args.output, 'w') as fh:
        fh.write(text)

    print(f"vespaw convert: wrote {args.output} "
          f"({len(atoms)} atoms, elements {', '.join(elements)})")
    for el, tag in used_tags:
        print(f"  {el:>2} -> !SPECIES NAME='{tag}'")


if __name__ == '__main__':
    main()
