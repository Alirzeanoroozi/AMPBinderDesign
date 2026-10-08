"""Shared ranking helpers for the structure panel.

Multimetric evaluation of a SHORT PEPTIDE bound to a large receptor. No single
confidence score is treated as a predicted affinity: PAE, pLDDT, ipSAE and
ipTM are confidence/geometry measures, not Kd.

Primary evidence (stage 4b, derived from the CIF/PAE Boltz-2 already wrote):
  * ipae_tgt_aligned  - interface PAE with the TARGET aligned and the PEPTIDE
    scored. Direction verified from ipsae.py's own byres header
    (i, AlignChn, ScoredChain): PAE row = aligned residue, column = scored.
    Unlike ipSAE in this direction it has no d0 floor, so it discriminates
    for 12-45 aa binders.
  * binder_iplddt - pLDDT of the peptide's interface residues. Complex pLDDT
    is useless here: it is dominated by the ~230 aa receptor (median 93-95
    regardless of how badly the peptide is placed).
  * catalytic/epitope contacts - binding the right site.
  * contact_density / buried_frac_binder - orthogonal interface geometry from
    coordinates alone, independent of any confidence head.

Secondary: ipTM (complex confidence only) and ipSAE_pep_aligned (the standard
ipSAE pair summary; demoted because it is the peptide-aligned/target-scored
direction and so can conceal uncertainty in PEPTIDE placement).

Deliberately NOT ranking terms: ipSAE_d0chn (dominated by full receptor
length), ipSAE_d0dom (no established peptide-specific standard), and
ipSAE_min/ipSAE_max (unreliable direction labels - across this run the
target-aligned direction is the minimum only ~59% of the time, ~30% ties).
All are retained as reported diagnostics.

Pose convergence across models/seeds is NOT included: only model_0 exists for
every prediction (single diffusion sample), so it cannot be computed without
new GPU inference.

Percentile ranks are computed WITHIN (target, peptide-length bin) because the
metrics are strongly length-confounded in opposite directions: median
ipae_tgt_aligned worsens 7.2 -> 14.2 A and median ipSAE_tgt_aligned collapses
0.030 -> 0.000 from 12-17 aa to 36-45 aa, while median binder pLDDT improves
42.5 -> 65.7. Unstratified ranking would systematically favour one length
class. Set STRATIFY_BY_LENGTH = False to disable.
"""
from __future__ import annotations

from difflib import SequenceMatcher

import pandas as pd

TARGETS = ("NDM5", "KPC3")
COLORS = {"NDM5": "#1f77b4", "KPC3": "#ff7f0e"}

# Percentile-rank weights (within target). Higher weight = more influence.
# Ranking is interface quality + epitope coverage only, plus developability
# penalties. Cationic/amphipathic "delivery proxy" terms (delivery_proxy,
# hydrophobic_moment, gravy) were deliberately removed: they are the biophysics
# of a membrane-lytic AMP, so weighting them reimposes the AMP objective on an
# enzyme-inhibitor design task. They remain in NUMERIC_COLS as annotations.
RANK_HIGHER = (
    ("binder_iplddt", 2.0),          # primary: peptide interface confidence
    ("n_catalytic_contacts", 1.5),   # primary: right site
    ("epitope_recall", 1.0),
    ("contact_density", 1.0),        # orthogonal geometry
    ("boltz2_iptm", 0.75),           # secondary complex confidence
    ("ipSAE_pep_aligned", 0.5),      # standard ipSAE summary, secondary
    ("buried_frac_binder", 0.5),
    ("interface_precision", 0.5),
)
RANK_LOWER = (
    ("ipae_tgt_aligned", 2.0),       # primary: peptide placement error, in A
    ("macrel_hemo_prob", 0.8),
    ("aggregation_proxy", 0.4),
    ("n_liabilities", 0.3),
)

# Peptide-length bins for stratified percentile ranking.
STRATIFY_BY_LENGTH = True
LENGTH_BINS = (11, 17, 23, 29, 35, 45)

# Structure gates. Thresholds are established conventions, not tuned:
#   50  = AlphaFold's "very low confidence" pLDDT boundary
#   0.5 = at least half the peptide engaged with the receptor
#
# ipae_tgt_aligned is LENGTH-DEPENDENT, so its threshold must be re-derived
# whenever BINDER_LEN changes. The original 10 A came from ipSAE's own PAE
# cutoff and was calibrated on the 12-45 aa pool, where it sat at the 41.1st
# percentile. At 50-100 aa the same percentile is 18.5 A (medians move
# 11.0 -> 18.1 A for NDM5 and 11.9 -> 20.7 A for KPC3), so a 10 A cut admitted
# 5/51 NDM5 and 0/45 KPC3 designs - it had silently become a near-total
# exclusion rather than a quality filter.
#
# 18.5 A is therefore the length-regime-matched equivalent of the original cut,
# NOT a loosening to make the funnel look better. Re-derive it the same way
# (match the percentile, don't guess) if the length range changes again.
DEFAULT_MIN_BINDER_IPLDDT = 50.0
DEFAULT_MAX_IPAE_TGT = 18.5
DEFAULT_MIN_BURIED_FRAC = 0.5
TOXIN_SAFE_WEIGHT = 1.0

DEFAULT_MIN_IPTM = 0.5
DEFAULT_N = 25
DEFAULT_MAX_IDENTITY = 0.8

NUMERIC_COLS = [c for c, _ in RANK_HIGHER] + [c for c, _ in RANK_LOWER] + [
    "boltz2_ptm",
    "boltz2_plddt",
    "ipSAE",
    "ipSAE_min",
    "ipSAE_max",
    "ipSAE_tgt_aligned",
    "ipSAE_d0chn",
    "ipSAE_d0dom",
    "pDockQ",
    "pDockQ2",
    "ipae_pep_aligned",
    "ipae_tgt_aligned_best",
    "binder_plddt_mean",
    "binder_plddt_min",
    "target_plddt_mean",
    "n_contacts",
    "n_if_pairs",
    "LIS",
    "length",
    "net_charge_pH7.4",
    "ampscanner_prob",
    "macrel_amp_prob",
    "hydramp_amp_prob",
    "hydramp_mic_prob",
    "toxinpred_hybrid_score",
    "delivery_proxy",
    "hydrophobic_moment",
    "gravy",
]


def find_col(columns, *names: str) -> str | None:
    cols = list(columns)
    for name in names:
        for c in cols:
            if c == name or c.endswith("." + name):
                return c
    return None


def is_true(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    return str(v).strip().lower() in ("true", "1", "yes")


def as_float(v, default=None):
    try:
        if v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def identity(a: str, b: str) -> float:
    return SequenceMatcher(None, a, b).ratio()


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in NUMERIC_COLS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def add_rank_score(df: pd.DataFrame) -> pd.DataFrame:
    """Add rank_score and binder_rank (1 = best) within each target."""
    out = coerce_numeric(df)
    out["toxin_safe"] = out.get("toxinpred_class", pd.Series("", index=out.index)).eq("Non-Toxin").astype(float)
    if "catalytic_ok" in out.columns:
        out["catalytic_ok_bool"] = out["catalytic_ok"].map(is_true)
    else:
        out["catalytic_ok_bool"] = False

    parts = []
    weights = []
    if STRATIFY_BY_LENGTH and "binder_len" in out.columns:
        out["_len_bin"] = pd.cut(
            pd.to_numeric(out["binder_len"], errors="coerce"),
            bins=list(LENGTH_BINS), labels=False, include_lowest=True,
        ).fillna(-1).astype(int)
        keys = ["target", "_len_bin"]
    else:
        out["_len_bin"] = -1
        keys = ["target"]
    grouped = out.groupby(keys, group_keys=False)

    for col, w in RANK_HIGHER:
        if col not in out.columns:
            continue
        score = grouped[col].transform(lambda s: s.rank(pct=True, na_option="keep")).fillna(0.0)
        parts.append(score)
        weights.append(w)
    for col, w in RANK_LOWER:
        if col not in out.columns:
            continue
        score = grouped[col].transform(
            lambda s: s.rank(pct=True, ascending=False, na_option="keep")
        ).fillna(0.0)
        parts.append(score)
        weights.append(w)
    if TOXIN_SAFE_WEIGHT and "toxin_safe" in out.columns:
        parts.append(out["toxin_safe"])
        weights.append(TOXIN_SAFE_WEIGHT)

    if not parts:
        out["rank_score"] = 0.0
    else:
        w = pd.Series(weights, dtype=float)
        mat = pd.concat(parts, axis=1)
        out["rank_score"] = mat.mul(w.values, axis=1).sum(axis=1) / w.sum()

    out["binder_rank"] = out.groupby("target")["rank_score"].rank(ascending=False, method="first")
    return out


def structure_gates(
    df: pd.DataFrame,
    min_iptm: float = DEFAULT_MIN_IPTM,
    require_catalytic: bool = True,
    require_non_toxin: bool = True,
    min_binder_iplddt: float | None = DEFAULT_MIN_BINDER_IPLDDT,
    max_ipae_tgt: float | None = DEFAULT_MAX_IPAE_TGT,
    min_buried_frac: float | None = DEFAULT_MIN_BURIED_FRAC,
) -> pd.Series:
    gates = pd.Series(True, index=df.index)
    if require_catalytic:
        if "catalytic_ok_bool" in df.columns:
            gates &= df["catalytic_ok_bool"].fillna(False)
        elif "catalytic_ok" in df.columns:
            gates &= df["catalytic_ok"].map(is_true)
    if min_iptm is not None and "boltz2_iptm" in df.columns:
        gates &= df["boltz2_iptm"].fillna(0.0) >= min_iptm
    # Peptide-placement gates. A design can post a high ipTM/ipSAE while the
    # peptide itself is unresolved, so require the peptide interface to be at
    # least marginally confident, placed, and actually buried.
    if min_binder_iplddt is not None and "binder_iplddt" in df.columns:
        gates &= df["binder_iplddt"].fillna(0.0) >= min_binder_iplddt
    if max_ipae_tgt is not None and "ipae_tgt_aligned" in df.columns:
        gates &= df["ipae_tgt_aligned"].fillna(1e9) <= max_ipae_tgt
    if min_buried_frac is not None and "buried_frac_binder" in df.columns:
        gates &= df["buried_frac_binder"].fillna(0.0) >= min_buried_frac
    if require_non_toxin and "toxinpred_class" in df.columns:
        gates &= df["toxinpred_class"].eq("Non-Toxin")
    return gates


def select_diverse(
    df: pd.DataFrame,
    n: int = DEFAULT_N,
    max_identity: float = DEFAULT_MAX_IDENTITY,
) -> pd.DataFrame:
    """Greedy take of top rank_score rows with pairwise sequence identity < max_identity."""
    ordered = df.sort_values(["rank_score", "boltz2_iptm"], ascending=False)
    picked_idx: list = []
    seqs: list[str] = []
    for idx, rec in ordered.iterrows():
        seq = str(rec.get("sequence") or "").upper()
        if not seq:
            continue
        if any(identity(seq, s) >= max_identity for s in seqs):
            continue
        picked_idx.append(idx)
        seqs.append(seq)
        if len(picked_idx) >= n:
            break
    if not picked_idx:
        return df.iloc[0:0].copy()
    return ordered.loc[picked_idx].reset_index(drop=True)


def select_panel(
    df: pd.DataFrame,
    n: int = DEFAULT_N,
    max_identity: float = DEFAULT_MAX_IDENTITY,
    min_iptm: float = DEFAULT_MIN_IPTM,
    require_catalytic: bool = True,
    require_non_toxin: bool = True,
    min_binder_iplddt: float | None = DEFAULT_MIN_BINDER_IPLDDT,
    max_ipae_tgt: float | None = DEFAULT_MAX_IPAE_TGT,
    min_buried_frac: float | None = DEFAULT_MIN_BURIED_FRAC,
) -> pd.DataFrame:
    """Prefer gated rows; if fewer than n survive diversity, relax iPTM then toxin.

    catalytic_ok is never dropped unless require_catalytic is False.
    """
    scored = add_rank_score(df)
    pep = dict(min_binder_iplddt=min_binder_iplddt, max_ipae_tgt=max_ipae_tgt,
               min_buried_frac=min_buried_frac)
    scored["passes_structure_gates"] = structure_gates(
        scored, min_iptm=min_iptm, require_catalytic=require_catalytic,
        require_non_toxin=require_non_toxin, **pep
    )

    pools = [scored[scored["passes_structure_gates"]]]
    if min_iptm is not None:
        pools.append(
            scored[
                structure_gates(
                    scored, min_iptm=None, require_catalytic=require_catalytic,
                    require_non_toxin=require_non_toxin, **pep
                )
            ]
        )
    if require_non_toxin:
        pools.append(
            scored[structure_gates(scored, min_iptm=None, require_catalytic=require_catalytic,
                                   require_non_toxin=False, **pep)]
        )
    if not require_catalytic:
        pools.append(scored)

    selected_parts = []
    used_ids: set[str] = set()
    remaining = n
    for pool in pools:
        if remaining <= 0:
            break
        pool = pool[~pool["design_id"].astype(str).isin(used_ids)]
        extra = select_diverse(pool, n=remaining, max_identity=max_identity)
        if extra.empty:
            continue
        selected_parts.append(extra)
        used_ids.update(extra["design_id"].astype(str))
        remaining = n - len(used_ids)

    if not selected_parts:
        return scored.iloc[0:0].copy()
    out = pd.concat(selected_parts, ignore_index=True)
    out["panel_rank"] = range(1, len(out) + 1)
    return out
