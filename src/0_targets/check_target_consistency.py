#!/usr/bin/env python3
"""Fail loudly if a target's structure and design sequence disagree.

The two design arms read the target from different files:
    BoltzGen      targets/<T>/design_domain.fasta      (sequence)
    RFdiffusion   targets/<T>/structures/<T>_target.pdb (coordinates)

If those disagree, the arms design against different proteins and their
designs are not comparable. That is exactly what happened: the stage-0
scaffolds are NDM-1 (5YPM) and KPC-2 (3DW0), the FASTA carried the
NDM-5 / KPC-3 mutations, and the structure did not - at positions that sit
on the designed interface.

Exit codes: 0 consistent, 1 mismatch, 2 could not check.

Usage:
  python src/0_targets/check_target_consistency.py --targets NDM5 KPC3
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "ASN2": "N", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "PHE": "F",
}


def pdb_chain_seq(pdb: Path, chain: str) -> dict:
    """{author residue number -> one-letter code} for CA atoms."""
    out = {}
    for line in pdb.open():
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            out[int(line[22:26])] = AA3.get(line[17:20].strip(), "X")
    return out


def numbering_map(tdir: Path):
    """[(precursor_pos, scaffold_resnum)] from numbering_map.csv.

    Scaffold (= PDB author) numbering is NOT always the precursor numbering:
    for KPC3 they differ by +2 in the C-terminal region, which is why a naive
    positional alignment of design_domain.fasta against the structure fails
    there. Comparing through this map is the only correct route.
    """
    nm = tdir / "numbering_map.csv"
    if not nm.is_file():
        return None
    out = []
    with nm.open(newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                out.append((int(row["precursor_pos"]),
                            int(float(row["scaffold_resnum"]))))
            except (KeyError, TypeError, ValueError):
                continue
    return out or None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", nargs="+", default=["NDM5", "KPC3"])
    ap.add_argument("--chain", default="A")
    ap.add_argument("--warn-only", action="store_true",
                    help="report and exit 0 even on mismatch")
    args = ap.parse_args()

    worst = 0
    for t in args.targets:
        tdir = REPO / "targets" / t
        fa = tdir / "design_domain.fasta"
        # prefer an explicitly mutated structure if one exists
        cands = [tdir / "structures" / f"{t}_target_mutated.pdb",
                 tdir / "structures" / f"{t}_target.pdb"]
        pdb = next((c for c in cands if c.is_file()), None)
        if not fa.is_file() or pdb is None:
            print(f"[{t}] CANNOT CHECK: missing {fa if not fa.is_file() else cands[-1]}")
            worst = max(worst, 2)
            continue
        pre_fa = tdir / "precursor.fasta"
        nmap = numbering_map(tdir)
        struct = pdb_chain_seq(pdb, args.chain)
        tag = f"[{t}] {pdb.name} chain {args.chain}"
        if nmap is None or not pre_fa.is_file():
            print(f"{tag}: CANNOT CHECK: need numbering_map.csv and precursor.fasta")
            worst = max(worst, 2)
            continue
        pre = "".join(l.strip() for l in pre_fa.read_text().splitlines()
                      if not l.startswith(">"))
        diffs, absent = [], 0
        for pp, sc in nmap:
            if not 1 <= pp <= len(pre):
                continue
            got = struct.get(sc)
            if got is None:
                absent += 1
                continue
            if got != pre[pp - 1]:
                diffs.append((pp, sc, pre[pp - 1], got))
        if not diffs:
            print(f"{tag}: consistent with the design sequence "
                  f"({len(nmap) - absent} residues compared"
                  + (f", {absent} unresolved in structure" if absent else "") + ")")
            continue
        print(f"{tag}: {len(diffs)} MISMATCH(ES) vs the design sequence")
        for pp, sc, want, got in diffs:
            print(f"    precursor {pp} / author {sc}: structure={got}  design={want}")
        print(f"    -> the two design arms would target different proteins.")
        print(f"    -> fix: python src/0_targets/apply_mutations.py --target {t} "
              f"--accept-truncated --out targets/{t}/structures/{t}_target_mutated.pdb")
        worst = max(worst, 1)

    if worst and args.warn_only:
        print("\n(--warn-only: continuing despite mismatch)")
        return 0
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
