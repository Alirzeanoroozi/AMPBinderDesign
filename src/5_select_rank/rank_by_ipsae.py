#!/usr/bin/env python3
"""Rank designs by ipSAE alone, highest to lowest.

Priority order, as requested:
    1. ipSAE_tgt_aligned   (target aligned, peptide scored)  - most important
    2. ipSAE_pep_aligned   (peptide aligned, target scored)
    3. ipSAE_min
    4. ipSAE_max

Two orderings are produced so you can use either:

  ipsae_rank_lex       strict priority order: sort on ipSAE_tgt_aligned
                       descending, break ties on ipSAE_pep_aligned, then
                       ipSAE_min, then ipSAE_max. The output file is sorted
                       this way.
  ipsae_rank_composite weighted mean of the four percentile ranks
                       (0.55 / 0.25 / 0.10 / 0.10) so the lower-priority
                       terms still influence the order rather than acting
                       only as tie-breaks, which with float values they
                       almost never do.

CAVEATS - this is a single-metric view, provided for comparison with the
multimetric panel. Read it knowing:

  * ipSAE's d0 is floored at ~1.04 A for any residue count <= 27, and
    ipSAE_tgt_aligned derives d0 from the PEPTIDE length. For 12-45 aa
    binders that direction is pinned at the floor, so this ranking is
    strongly length-biased: median ipSAE_tgt_aligned falls 0.030 -> 0.000
    from the 12-17 aa bin to the 36-45 aa bin.
  * ipSAE_min / ipSAE_max are NOT direction labels. Across this run the
    target-aligned direction is the minimum only 59.4% of the time, with
    29.8% ties. They are included because they were requested, but
    ipSAE_tgt_aligned / ipSAE_pep_aligned are the meaningful pair.
  * ipSAE is a confidence score, not a predicted affinity or Kd.

Context columns (ipae_tgt_aligned, binder_iplddt, ...) are carried so a high
ipSAE design can be checked against whether the peptide is actually resolved
and placed.

Usage:
  python src/5_select_rank/rank_by_ipsae.py
  python src/5_select_rank/rank_by_ipsae.py --top 25
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TARGETS = ("NDM5", "KPC3")

IPSAE_PRIORITY = ["ipSAE_tgt_aligned", "ipSAE_pep_aligned", "ipSAE_min", "ipSAE_max"]
COMPOSITE_WEIGHTS = {"ipSAE_tgt_aligned": 0.55, "ipSAE_pep_aligned": 0.25,
                     "ipSAE_min": 0.10, "ipSAE_max": 0.10}

OUT_COLS = [
    "ipsae_rank_lex", "ipsae_rank_composite", "ipsae_composite",
    "design_id", "target", "source", "sequence", "length",
    # the requested ipSAE family, in priority order
    "ipSAE_tgt_aligned", "ipSAE_pep_aligned", "ipSAE_min", "ipSAE_max",
    "ipSAE", "ipSAE_d0chn", "ipSAE_d0dom",
    # context: is the peptide actually resolved and placed?
    "ipae_tgt_aligned", "binder_iplddt", "binder_plddt_mean",
    "boltz2_iptm", "epitope_recall", "n_catalytic_contacts", "catalytic_ok",
    "contact_density", "buried_frac_binder",
    # cross-reference against the multimetric panel
    "passes_structure_gates", "in_multimetric_panel",
    "multimetric_rank_score", "multimetric_binder_rank",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default=str(REPO / "results"))
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--top", type=int, default=0,
                    help="also write a top-N file (0 = full ranking only)")
    ap.add_argument("--gated-only", action="store_true",
                    help="restrict to designs passing the structure gates")
    args = ap.parse_args()
    R = Path(args.results_dir)

    for t in args.targets:
        src = R / f"ranked_{t}.csv"
        if not src.is_file():
            print(f"[{t}] missing {src} - run stage 5 first")
            continue
        df = pd.read_csv(src, low_memory=False)

        panel = set()
        sel = R / f"selected_{t}.csv"
        if sel.is_file():
            panel = set(pd.read_csv(sel, low_memory=False)["design_id"].astype(str))
        df["in_multimetric_panel"] = df["design_id"].astype(str).isin(panel)
        df = df.rename(columns={"rank_score": "multimetric_rank_score",
                                "binder_rank": "multimetric_binder_rank"})

        for c in IPSAE_PRIORITY:
            if c not in df.columns:
                df[c] = pd.NA
            df[c] = pd.to_numeric(df[c], errors="coerce")

        before = len(df)
        df = df[df["ipSAE_tgt_aligned"].notna()].copy()
        dropped = before - len(df)
        if args.gated_only and "passes_structure_gates" in df.columns:
            df = df[df["passes_structure_gates"].astype(str) == "True"].copy()

        # 1. strict priority order
        df = df.sort_values(IPSAE_PRIORITY, ascending=False, kind="mergesort")
        df["ipsae_rank_lex"] = range(1, len(df) + 1)

        # 2. weighted composite of percentile ranks
        comp = sum(df[c].rank(pct=True) * w for c, w in COMPOSITE_WEIGHTS.items())
        df["ipsae_composite"] = comp.round(6)
        df["ipsae_rank_composite"] = (
            df["ipsae_composite"].rank(ascending=False, method="first").astype(int)
        )

        cols = [c for c in OUT_COLS if c in df.columns]
        out = R / f"ipsae_ranked_{t}.csv"
        df[cols].to_csv(out, index=False)
        print(f"[{t}] {len(df)} designs ranked by ipSAE "
              f"({dropped} dropped: no ipSAE value) -> {out.relative_to(REPO)}")
        if args.top:
            top = R / f"ipsae_top{args.top}_{t}.csv"
            df[cols].head(args.top).to_csv(top, index=False)
            print(f"      top {args.top} -> {top.relative_to(REPO)}")
        agree = df.head(25)["in_multimetric_panel"].sum()
        print(f"      of the ipSAE top-25, {agree} are also in the multimetric panel")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
