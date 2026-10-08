#!/usr/bin/env python3
"""Apply the documented target mutations to the design structure.

Why this exists
---------------
Stage 0 builds each target from a crystallographic scaffold that is NOT the
variant we are designing against:

    NDM5  scaffold 5YPM chain A  ->  needs V88L, M154L   (NDM-1 -> NDM-5)
    KPC3  scaffold 3DW0 chain A  ->  needs H274Y         (KPC-2 -> KPC-3)

`design_domain.fasta` already carries the mutated (correct) sequence, and
`prep_report.md` says "apply to structure with apply_mutations.py if needed" -
but that script was never written, so the mutations were never applied to the
coordinates. BoltzGen reads the FASTA so it designed against the right
sequence; RFdiffusion reads the PDB, so it designed against the WRONG residue
at those positions. All three sites sit at the designed interface:

    NDM5 V88   3.6 A from the nearest epitope residue
    NDM5 M154  3.7 A
    KPC3 H274  0.0 A - it IS an epitope residue

Sidechain handling
------------------
No sidechain rebuilding library is available in these environments, and this
script will not invent coordinates. It renames the residue and TRUNCATES the
sidechain to backbone + CB, which is always geometrically valid for the new
residue type. Atoms beyond CB are removed and must be rebuilt by a packer
(PyMOL mutagenesis, Rosetta, pdbfixer, FoldX) before the structure is used
for design. The script prints exactly what it removed and exits non-zero
unless --accept-truncated is passed, so a run cannot silently proceed on a
half-built sidechain.

Usage:
  python src/0_targets/apply_mutations.py --target NDM5 --dry-run
  python src/0_targets/apply_mutations.py --target NDM5 --accept-truncated \
      --out targets/NDM5/structures/NDM5_target_mutated.pdb
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
AA1 = {v: k for k, v in AA3.items()}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
KEEP = BACKBONE | {"CB"}            # CB is valid for every residue except GLY
MUT_RE = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])\b")


def mutations_from_report(report: Path) -> list[tuple[str, int, str]]:
    """Parse 'Mutations vs template ...: V88L, M154L' out of prep_report.md."""
    if not report.is_file():
        return []
    for line in report.read_text().splitlines():
        if "mutation" in line.lower() and "template" in line.lower():
            return [(m.group(1), int(m.group(2)), m.group(3))
                    for m in MUT_RE.finditer(line.split(":", 1)[-1])]
    return []


def _rel(p: Path) -> str:
    """Display path relative to REPO when possible, else as given."""
    try:
        return str(Path(p).resolve().relative_to(REPO))
    except ValueError:
        return str(p)


def chain_seq(pdb: Path, chain: str) -> dict[int, str]:
    out = {}
    for line in pdb.open():
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            out[int(line[22:26])] = AA3.get(line[17:20].strip(), "X")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--pdb", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--mutations", nargs="*", default=None,
                    help="override, e.g. V88L M154L (default: read prep_report.md)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--accept-truncated", action="store_true",
                    help="acknowledge that sidechains beyond CB are removed and "
                         "must be rebuilt before use")
    args = ap.parse_args()

    tdir = REPO / "targets" / args.target
    pdb = Path(args.pdb) if args.pdb else tdir / "structures" / f"{args.target}_target.pdb"
    if not pdb.is_file():
        sys.exit(f"no such structure: {pdb}")

    if args.mutations:
        muts = [(m[0], int(m[1:-1]), m[-1]) for m in args.mutations]
    else:
        muts = mutations_from_report(tdir / "prep_report.md")
    if not muts:
        print(f"[{args.target}] no mutations documented - nothing to do")
        return 0

    seq = chain_seq(pdb, args.chain)
    todo, already, bad = [], [], []
    for frm, pos, to in muts:
        have = seq.get(pos)
        if have is None:
            bad.append((frm, pos, to, "position absent from structure"))
        elif have == to:
            already.append((frm, pos, to))
        elif have != frm:
            bad.append((frm, pos, to, f"structure has {have}, expected {frm}"))
        else:
            todo.append((frm, pos, to))

    print(f"[{args.target}] {_rel(pdb)} chain {args.chain}")
    for frm, pos, to in already:
        print(f"  {frm}{pos}{to}: already applied")
    for frm, pos, to, why in bad:
        print(f"  {frm}{pos}{to}: CANNOT APPLY - {why}")
    for frm, pos, to in todo:
        print(f"  {frm}{pos}{to}: will mutate {AA1[frm]} -> {AA1[to]}, "
              f"truncating sidechain to backbone+CB")
    if bad:
        return 2
    if not todo:
        print("  nothing to change")
        return 0
    if args.dry_run:
        print("  (dry run - no file written)")
        return 0
    if not args.accept_truncated:
        print("\n  REFUSING to write: the new sidechains would be incomplete "
              "(backbone+CB only).\n  Rebuild them with a packer, or pass "
              "--accept-truncated if you understand this.", file=sys.stderr)
        return 3

    out = Path(args.out).resolve() if args.out else pdb.with_name(f"{pdb.stem}_mutated.pdb")
    bypos = {p: (f, t) for f, p, t in todo}
    removed, kept_lines = [], []
    for line in pdb.open():
        if line.startswith(("ATOM", "HETATM")) and line[21] == args.chain:
            pos = int(line[22:26])
            if pos in bypos:
                atom = line[12:16].strip()
                if atom not in KEEP:
                    removed.append((pos, atom))
                    continue
                line = line[:17] + AA1[bypos[pos][1]].ljust(3) + line[20:]
        kept_lines.append(line)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(kept_lines))
    print(f"\n  wrote {_rel(out)}")
    print(f"  removed {len(removed)} sidechain atoms: "
          + ", ".join(f"{p}:{a}" for p, a in removed))
    print("  THESE SIDECHAINS ARE INCOMPLETE - rebuild before using for design.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
