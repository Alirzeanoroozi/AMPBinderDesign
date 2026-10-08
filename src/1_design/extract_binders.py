#!/usr/bin/env python3
"""Extract binder sequences from generated complexes.

Only *finished designs* become candidates. A design-output tree holds several
overlapping copies of the same run, and most of them are not designs:

  boltzgen_outputs/<T>/intermediate_designs/                 backbone placeholders,
                                                             NOT sequence-designed
  boltzgen_outputs/<T>/intermediate_designs_inverse_folded/  real designs
  .../intermediate_designs_inverse_folded/refold_cif/        same designs again
  .../intermediate_designs_inverse_folded/refold_design_cif/ same designs again
  boltzgen_outputs/<T>/final_ranked_designs/final_200_designs/   BoltzGen's ranked output
  .../final_200_designs/before_refolding/                    same designs again
  .../intermediate_ranked_10_designs/                        subset of the above
  boltzgen_outputs/<T>/<T>.cif                               the input template
  rfdiffusion_outputs/<T>_inverse_folding/                   real designs

Recursing the whole tree (the old behaviour) swept all of it into one pool, so
pre-inverse-folding placeholders competed as candidates and every real design
was entered 3-4 times under different ids. This module classifies each file
explicitly and keeps one record per distinct binder sequence.

BoltzGen's own QC
-----------------
`boltzgen run` does not stop at inverse folding: its filter step refolds every
design (in complex, design-chain-only, and in isolation), applies hard
thresholds, then ranks by quality and diversity down to --budget. It rejects
89.5% (NDM5) and 96.0% (KPC3) of its own inverse-folding output. That verdict
is recorded per design in
    boltzgen_outputs/<T>/final_ranked_designs/all_designs_metrics.csv
as `pass_filters`. Harvesting intermediate_designs_inverse_folded wholesale
therefore imports the ~90-96% BoltzGen had already judged unusable, which is
where the low-complexity sequences came from.

So boltzgen_ifold designs are admitted ONLY if BoltzGen marked them
pass_filters=True. That yields the full QC-passing set (525 NDM5 / 201 KPC3)
rather than just the diversity-selected top --budget, while excluding the
designs BoltzGen rejected.

The RFdiffusion arm originally had NO quality control of any kind, so it would
have outnumbered the QC-passed BoltzGen arm ~20-50x with unselected designs.
scripts/run_1d_rfdiffusion_qc.slurm now runs the SAME boltzgen filter steps
(folding, design_folding, analysis, filtering) over the RFdiffusion
inverse-folding output, and rfdiffusion_ifold is gated on that verdict too.
Both arms therefore face identical criteria, including the refold-in-isolation
check that nothing downstream can replicate. Pass --no-boltzgen-qc to disable
gating for both.

Binder = shortest protein chain unless --binder-chain is supplied.

Usage (from AMPBinderDesign):
  conda activate ampbinder
  python src/1_design/extract_binders.py
"""
from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "KCX": "K",
    "UNK": "X",
}

# Preference order when the same sequence is produced by more than one source.
SOURCE_PRIORITY = ("boltzgen_final", "boltzgen_ifold", "rfdiffusion_ifold")

# Any binder longer than this is assumed to be a misidentified target chain.
DEFAULT_MAX_BINDER_LEN = 120


# Where each QC-gated source's BoltzGen filter verdict lives, relative to REPO.
# Both arms are gated on the SAME criteria: the RFdiffusion designs are put
# through boltzgen's folding / design_folding / analysis / filtering steps by
# scripts/run_1d_rfdiffusion_qc.slurm, which writes an all_designs_metrics.csv
# in exactly the same format as the BoltzGen arm's.
QC_METRICS = {
    "boltzgen_final": "boltzgen_outputs/{target}/final_ranked_designs/all_designs_metrics.csv",
    "boltzgen_ifold": "boltzgen_outputs/{target}/final_ranked_designs/all_designs_metrics.csv",
    "rfdiffusion_ifold": "rfdiffusion_qc/{target}/final_ranked_designs/all_designs_metrics.csv",
}

# boltzgen_final was originally treated as QC-exempt on the assumption that
# final_<budget>_designs/ only ever contains designs passing the filters. That
# is FALSE: boltzgen writes `budget` designs regardless of how many passed, so
# whenever fewer designs pass than the budget the directory is padded with
# REJECTED designs. Measured on a 200-design run with budget 20: only 6 (NDM5)
# and 8 (KPC3) passed, yet 20 were written - so 14/20 and 12/20 had failed QC,
# including 5 and 2 that carried a cysteine. It is now gated like every other
# source.
#
# Its filenames carry a rank prefix (rank001_NDM5_026) while the metrics table
# keys on the bare design id (NDM5_026), so the prefix is stripped FOR THE QC
# LOOKUP ONLY. The design_id itself keeps the prefix - stripping it there is the
# stage-4 bug that silently discarded every ranked design's scores.
_QC_ID_STRIP = {"boltzgen_final": re.compile(r"^rank\d+_", re.IGNORECASE)}


def qc_lookup_id(source: str, stem: str) -> str:
    rx = _QC_ID_STRIP.get(source)
    return rx.sub("", stem) if rx else stem


def qc_pass_set(path: str) -> Optional[set]:
    """design ids a boltzgen filter step marked pass_filters=True, or None."""
    if not os.path.isfile(path):
        return None
    ok: set = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if str(row.get("pass_filters", "")).strip().lower() in ("true", "1", "1.0"):
                did = (row.get("id") or "").strip()
                if did:
                    ok.add(did)
    return ok


def chain_sequence(chain) -> str:
    letters = []
    for res in chain:
        if res.id[0] != " ":
            continue
        letters.append(THREE_TO_ONE.get(res.get_resname().strip(), "X"))
    return "".join(letters)


def load_structure(path: str, cif_parser: MMCIFParser, pdb_parser: PDBParser):
    if path.endswith(".cif"):
        return cif_parser.get_structure(os.path.basename(path), path)
    if path.endswith(".pdb"):
        return pdb_parser.get_structure(os.path.basename(path), path)
    raise ValueError(f"unsupported structure extension: {path}")


def binder_sequence(
    structure_path: str,
    cif_parser: MMCIFParser,
    pdb_parser: PDBParser,
    binder_chain: Optional[str] = None,
) -> Tuple[str, str]:
    """Return (binder_seq, chain_id) using the shortest amino-acid chain."""
    structure = load_structure(structure_path, cif_parser, pdb_parser)
    model = next(structure.get_models())
    chains: List[Tuple[str, str]] = []
    for chain in model:
        seq = chain_sequence(chain)
        if seq:
            chains.append((chain.id, seq))
    if not chains:
        raise ValueError(f"no protein chains in {structure_path}")
    if binder_chain is not None:
        for chain_id, seq in chains:
            if chain_id == binder_chain:
                return seq, chain_id
        raise ValueError(f"binder chain {binder_chain!r} not found in {structure_path}")
    chain_id, seq = min(chains, key=lambda x: len(x[1]))
    return seq, chain_id


def target_from_path(path: str, targets: List[str]) -> Optional[str]:
    """Assign a CIF to a target if NDM5 or KPC3 appears in the file path."""
    blob = path.replace("\\", "/").upper()
    hits = [t for t in targets if t.upper() in blob]
    if not hits:
        return None
    return max(hits, key=len)


def classify(rel_path: str, targets: List[str]) -> Tuple[Optional[str], str]:
    """Map a repo-relative path to (source_label, reject_reason).

    source_label is None when the file is not a finished design; reject_reason
    then says why. Exclusions are tested before inclusions so that a copy
    directory nested under a kept directory is still dropped.
    """
    rel = rel_path.replace("\\", "/")
    padded = "/" + rel.strip("/")
    base = padded.rsplit("/", 1)[-1]
    stem = base.rsplit(".", 1)[0]

    if "/traj/" in padded:
        return None, "trajectory"
    if "_native" in base:
        return None, "native_reference"
    if "/boltzgen_outputs/" in padded and stem in targets:
        return None, "input_template"
    if "/intermediate_designs/" in padded:
        return None, "placeholder_pre_inverse_folding"
    if "/before_refolding/" in padded:
        return None, "duplicate_before_refolding"
    if "/refold_cif/" in padded or "/refold_design_cif/" in padded:
        return None, "duplicate_refold_copy"
    if "/intermediate_ranked_10_designs/" in padded:
        return None, "duplicate_ranked_subset"

    if "/final_ranked_designs/" in padded:
        return "boltzgen_final", ""
    if "/intermediate_designs_inverse_folded/" in padded:
        return "boltzgen_ifold", ""
    if "_inverse_folding/" in padded:
        return "rfdiffusion_ifold", ""
    return None, "unclassified_directory"


def list_structures(directory: str) -> List[str]:
    """All .cif/.pdb files under directory, recursively."""
    out = []
    if not os.path.isdir(directory):
        return out
    for root, _dirs, files in os.walk(directory):
        for name in files:
            if name.endswith(".cif") or name.endswith(".pdb"):
                out.append(os.path.join(root, name))
    out.sort()
    return out


def extract_all(
    directories: List[str],
    targets: List[str],
    binder_chain: Optional[str] = None,
    sources: Optional[List[str]] = None,
    max_binder_len: int = DEFAULT_MAX_BINDER_LEN,
    allow_unclassified: bool = False,
    use_boltzgen_qc: bool = True,
) -> Tuple[Dict[str, List[SeqRecord]], Dict[str, List[dict]]]:
    cif_parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)
    wanted = set(sources or SOURCE_PRIORITY)

    # qc[(target, source)] -> set of ids passing that source's QC, or None
    qc: Dict[tuple, Optional[set]] = {}
    if use_boltzgen_qc:
        for t in targets:
            for src, tmpl in QC_METRICS.items():
                path = os.path.join(REPO, tmpl.format(target=t))
                got = qc_pass_set(path)
                qc[(t, src)] = got
                if got is None:
                    print(f"  [{t}] WARNING: {src} is NOT QC-gated "
                          f"(missing {os.path.relpath(path, REPO)})")
                else:
                    print(f"  [{t}] {src}: {len(got)} designs passed boltzgen QC")

    paths = []
    for directory in directories:
        paths.extend(list_structures(directory))
    paths = sorted(set(paths))
    print(f"Found {len(paths)} structures under {', '.join(directories)}")

    rejected: Counter = Counter()
    skipped = 0
    unmatched = 0
    # (target, sequence) -> provenance dict, best source wins
    best: Dict[Tuple[str, str], dict] = {}
    collisions: Dict[str, set] = {}

    for i, path in enumerate(paths, 1):
        rel = os.path.relpath(path, REPO)
        target = target_from_path(path, targets)
        if target is None:
            unmatched += 1
            continue
        source, reason = classify(rel, targets)
        if source is None:
            if reason == "unclassified_directory" and allow_unclassified:
                source = "unclassified"
            else:
                rejected[reason] += 1
                continue
        if source not in wanted and source != "unclassified":
            rejected[f"source_not_selected:{source}"] += 1
            continue
        # Honour the boltzgen filter verdict for every QC-gated source rather
        # than re-deriving a worse one downstream. Applied identically to both
        # design arms so neither enters the pool unselected.
        passed = qc.get((target, source))
        if passed is not None:
            stem_id = os.path.splitext(os.path.basename(path))[0]
            if qc_lookup_id(source, stem_id) not in passed:
                rejected[f"qc_rejected:{source}"] += 1
                continue

        try:
            seq, chain_id = binder_sequence(path, cif_parser, pdb_parser, binder_chain)
        except Exception as exc:
            print(f"  skip {rel}: {exc}")
            skipped += 1
            continue
        if not seq:
            print(f"  skip {rel}: empty binder sequence")
            skipped += 1
            continue
        if len(seq) > max_binder_len:
            rejected[f"binder_longer_than_{max_binder_len}"] += 1
            continue

        design_id = os.path.splitext(os.path.basename(path))[0]
        collisions.setdefault(f"{target}/{design_id}", set()).add(source)
        rec = {
            "design_id": design_id,
            "target": target,
            "source": source,
            "boltzgen_qc": ("pass" if qc.get((target, source)) is not None
                            else ("ungated" if source in QC_METRICS else "n/a")),
            "chain": chain_id,
            "length": len(seq),
            "sequence": seq,
            "file": rel,
        }
        key = (target, seq)
        prior = best.get(key)
        if prior is None or _rank(source) < _rank(prior["source"]) or (
            _rank(source) == _rank(prior["source"]) and design_id < prior["design_id"]
        ):
            best[key] = rec
        if i % 500 == 0 or i == len(paths):
            print(f"  parsed {i}/{len(paths)}")

    # A design_id must name exactly one design. Fail loudly rather than suffix.
    clashes = {k: v for k, v in collisions.items() if len(v) > 1}
    if clashes:
        raise SystemExit(
            "design_id collision across sources (the id would no longer be unique): "
            + ", ".join(f"{k} <- {sorted(v)}" for k, v in sorted(clashes.items())[:10])
        )

    records: Dict[str, List[SeqRecord]] = {t: [] for t in targets}
    provenance: Dict[str, List[dict]] = {t: [] for t in targets}
    for rec in sorted(best.values(), key=lambda r: (r["target"], _rank(r["source"]), r["design_id"])):
        records[rec["target"]].append(
            SeqRecord(
                Seq(rec["sequence"]),
                id=rec["design_id"],
                description=(
                    f"target={rec['target']} source={rec['source']} chain={rec['chain']} "
                    f"length={rec['length']} file={rec['file']}"
                ),
            )
        )
        provenance[rec["target"]].append(rec)

    kept = sum(len(v) for v in records.values())
    print(f"\nkept {kept} distinct binder sequences; skipped={skipped} unmatched={unmatched}")
    if rejected:
        print("rejected (not finished designs):")
        for reason, n in rejected.most_common():
            print(f"  {n:7d}  {reason}")
    for target in targets:
        by_src = Counter(r["source"] for r in provenance[target])
        print(f"  [{target}] {len(provenance[target])} designs: {dict(by_src)}")
    return records, provenance


def _rank(source: str) -> int:
    try:
        return SOURCE_PRIORITY.index(source)
    except ValueError:
        return len(SOURCE_PRIORITY)


def write_fasta(records: List[SeqRecord], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = SeqIO.write(records, path, "fasta")
    print(f"  wrote {n} sequences -> {path}")


def write_provenance(rows: List[dict], path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cols = ["design_id", "target", "source", "boltzgen_qc", "chain", "length",
            "sequence", "file"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"  wrote {len(rows)} rows -> {path}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--outputs-dir",
        nargs="+",
        default=[os.path.join(REPO, "boltzgen_outputs")],
        help="one or more roots to recurse for CIF/PDB files",
    )
    ap.add_argument(
        "--out-dir",
        default=os.path.join(REPO, "fastas"),
        help="directory for NDM5.fasta and KPC3.fasta",
    )
    ap.add_argument("--targets", nargs="+", default=["NDM5", "KPC3"])
    ap.add_argument(
        "--binder-chain",
        default=None,
        help="extract this chain instead of choosing the shortest protein chain",
    )
    ap.add_argument(
        "--sources",
        nargs="+",
        default=list(SOURCE_PRIORITY),
        choices=list(SOURCE_PRIORITY),
        help="which design sources to keep (default: all finished-design sources)",
    )
    ap.add_argument(
        "--max-binder-len",
        type=int,
        default=DEFAULT_MAX_BINDER_LEN,
        help="reject a 'binder' longer than this; guards against picking a target chain",
    )
    ap.add_argument(
        "--no-boltzgen-qc",
        action="store_true",
        help="do NOT gate boltzgen_ifold on BoltzGen's own pass_filters verdict "
             "(restores the old behaviour of importing its pre-QC output)",
    )
    ap.add_argument(
        "--allow-unclassified",
        action="store_true",
        help="also keep structures in directories this script does not recognise",
    )
    args = ap.parse_args()

    by_target, provenance = extract_all(
        args.outputs_dir,
        args.targets,
        args.binder_chain,
        sources=args.sources,
        max_binder_len=args.max_binder_len,
        allow_unclassified=args.allow_unclassified,
        use_boltzgen_qc=not args.no_boltzgen_qc,
    )
    for target, recs in by_target.items():
        write_fasta(recs, os.path.join(args.out_dir, f"{target}.fasta"))
        write_provenance(provenance[target], os.path.join(args.out_dir, f"{target}_provenance.csv"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
