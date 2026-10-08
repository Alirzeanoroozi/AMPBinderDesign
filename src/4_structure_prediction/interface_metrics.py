#!/usr/bin/env python3
"""Peptide-aware interface metrics from EXISTING Boltz-2 outputs. No GPU.

Post-processing only: reads the model CIF and pae_*.npz that Boltz-2 already
wrote, plus the per-complex ipsae.py .txt, and derives metrics that are
meaningful for a short peptide bound to a large receptor.

Why these and not ipSAE alone
-----------------------------
ipSAE's d0 is floored at ~1.04 A for any residue count <= 27, and the
target-aligned direction derives d0 from the PEPTIDE length, so for a 12-45 aa
binder that direction is pinned at the floor and cannot discriminate. The raw
interface PAE in the same direction has no such floor, so it is used directly.
Complex pLDDT is likewise uninformative here: it is dominated by the ~230 aa
receptor, so the peptide chain is scored separately.

Direction convention (verified, not assumed)
--------------------------------------------
ipsae.py's own by-residue header reads
    i, AlignChn, ScoredChain, AlignResNum, ...
and it indexes pae_matrix[i] for i in chain1, selecting columns in chain2.
Therefore PAE row = ALIGNED residue, PAE column = SCORED residue, and the
.txt row with Chn1 == <target chain> is "target aligned, peptide scored".
Both directions are emitted under explicit names; min/max are NOT used as
direction labels because across this run the target-aligned direction is the
minimum only ~59% of the time (~30% ties).

Nothing here is an affinity estimate. PAE, pLDDT, ipSAE and ipTM are
confidence and geometry measures, not predicted Kd.

Usage:
  python src/4_structure_prediction/interface_metrics.py --workers 16
"""
from __future__ import annotations

import argparse
import csv
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TARGETS = ("NDM5", "KPC3")
_MODEL_SUFFIX_RE = re.compile(r"_model_\d+$", re.IGNORECASE)
_SKIP_DIR_PARTS = {"processed", "mols", "msa"}

COLS = [
    "design_id",
    "target_chain", "binder_chain", "target_len", "binder_len",
    # primary: peptide placement confidence, no d0 floor
    "ipae_tgt_aligned", "ipae_pep_aligned", "ipae_tgt_aligned_best",
    # primary: peptide-chain pLDDT (complex pLDDT is receptor-dominated)
    "binder_plddt_mean", "binder_plddt_min", "binder_iplddt", "target_plddt_mean",
    # explicitly-labelled directional ipSAE (from ipsae.py asym rows)
    "ipSAE_tgt_aligned", "ipSAE_pep_aligned",
    # orthogonal interface geometry (coordinates only, no confidence head)
    "n_if_pairs", "n_contacts", "contact_density", "buried_frac_binder",
    "status",
]


def parse_cif(path: Path):
    """Per-residue records in CIF order (== PAE/pLDDT array order)."""
    res, order = {}, []
    for line in path.open():
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        if len(f) < 18:
            continue
        ch, seq, el, atom = f[9], int(f[6]), f[2], f[3]
        key = (ch, seq)
        if key not in res:
            res[key] = {"chain": ch, "ca": None, "cb": None, "heavy": [],
                        "plddt": float(f[17])}
            order.append(key)
        xyz = (float(f[10]), float(f[11]), float(f[12]))
        if atom == "CA":
            res[key]["ca"] = xyz
        elif atom == "CB":
            res[key]["cb"] = xyz
        if el != "H":
            res[key]["heavy"].append(xyz)
    return [res[k] for k in order]


def directional_ipsae(txt: Path, target_chain: str, binder_chain: str):
    """ipSAE for each asym row, keyed by its ALIGN chain (Chn1)."""
    out = {}
    if not txt.exists():
        return out
    hdr = None
    for line in txt.read_text().splitlines():
        line = line.strip()
        if line.startswith("Chn1,"):
            hdr = [x.strip() for x in line.split(",")]
            continue
        if not hdr:
            continue
        parts = [x.strip() for x in line.split(",")]
        if len(parts) != len(hdr):
            continue
        r = dict(zip(hdr, parts))
        if r.get("Type") != "asym":
            continue
        if {r.get("Chn1"), r.get("Chn2")} != {target_chain, binder_chain}:
            continue
        try:
            v = float(r["ipSAE"])
        except (KeyError, ValueError):
            continue
        # Chn1 is the align chain (verified from ipsae.py's byres header)
        out["ipSAE_tgt_aligned" if r["Chn1"] == target_chain
            else "ipSAE_pep_aligned"] = v
    return out


def score_one(did: str, cif_s: str, pae_s: str, txt_s: str,
              dist_cutoff: float, contact_cutoff: float) -> dict:
    row = {c: "" for c in COLS}
    row["design_id"] = did
    try:
        cif, pae_p = Path(cif_s), Path(pae_s)
        rs = parse_cif(cif)
        pae = np.load(pae_p)["pae"]
        if len(rs) != pae.shape[0]:
            row["status"] = f"len_mismatch cif={len(rs)} pae={pae.shape[0]}"
            return row
        chains = np.array([r["chain"] for r in rs])
        uniq, counts = np.unique(chains, return_counts=True)
        if len(uniq) != 2:
            row["status"] = f"expected_2_chains_got_{len(uniq)}"
            return row
        # target = the longer chain; binder = the shorter. Verified A/B in this
        # run but derived from length so a flipped YAML cannot silently invert.
        tgt_c = uniq[int(np.argmax(counts))]
        pep_c = uniq[int(np.argmin(counts))]
        tgt = np.where(chains == tgt_c)[0]
        pep = np.where(chains == pep_c)[0]

        cb = np.array([r["cb"] or r["ca"] for r in rs], dtype=float)
        sub = np.linalg.norm(cb[tgt][:, None, :] - cb[pep][None, :, :], axis=2)
        pairs = sub < dist_cutoff
        if not pairs.any():
            row["status"] = "no_interface"
            return row

        block_t = pae[np.ix_(tgt, pep)]   # row=target(aligned), col=peptide(scored)
        block_p = pae[np.ix_(pep, tgt)]   # row=peptide(aligned), col=target(scored)
        row["ipae_tgt_aligned"] = round(float(block_t[pairs].mean()), 3)
        row["ipae_pep_aligned"] = round(float(block_p[pairs.T].mean()), 3)
        row["ipae_tgt_aligned_best"] = round(float(block_t[pairs].min()), 3)

        plddt = np.array([r["plddt"] for r in rs], dtype=float)
        pep_if = pep[pairs.sum(axis=0) > 0]
        row["binder_plddt_mean"] = round(float(plddt[pep].mean()), 2)
        row["binder_plddt_min"] = round(float(plddt[pep].min()), 2)
        row["binder_iplddt"] = round(float(plddt[pep_if].mean()), 2)
        row["target_plddt_mean"] = round(float(plddt[tgt].mean()), 2)

        ncon = 0
        for jj, j in enumerate(pep):
            close = tgt[sub[:, jj] < dist_cutoff + 4.0]
            if close.size == 0:
                continue
            b = np.asarray(rs[j]["heavy"], dtype=float)
            for i in close:
                a = np.asarray(rs[i]["heavy"], dtype=float)
                if np.min(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)) < contact_cutoff:
                    ncon += 1
        row["n_if_pairs"] = int(pairs.sum())
        row["n_contacts"] = ncon
        row["contact_density"] = round(ncon / len(pep), 3)
        row["buried_frac_binder"] = round(len(pep_if) / len(pep), 3)

        row["target_chain"] = tgt_c
        row["binder_chain"] = pep_c
        row["target_len"] = len(tgt)
        row["binder_len"] = len(pep)
        row.update(directional_ipsae(Path(txt_s), tgt_c, pep_c))
        row["status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        row["status"] = f"ERROR {type(exc).__name__}: {exc}"
    return row


def discover(boltz_out: Path, target: str):
    jobs = []
    for cif in sorted(boltz_out.rglob("*_model_0.cif")):
        if any(p in _SKIP_DIR_PARTS for p in cif.parts):
            continue
        did = _MODEL_SUFFIX_RE.sub("", cif.stem)
        if target and target not in did:
            continue
        pae = cif.with_name(f"pae_{cif.stem}.npz")
        if not pae.exists():
            continue
        jobs.append((did, str(cif), str(pae),
                     str(cif.with_name(f"{cif.stem}_10_10.txt"))))
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--boltz-out", default=str(REPO / "boltz_results"))
    ap.add_argument("--targets", nargs="+", default=list(TARGETS))
    ap.add_argument("--out-dir", default=str(REPO / "results"))
    ap.add_argument("--dist-cutoff", type=float, default=10.0,
                    help="CB-CB A defining an interface residue pair")
    ap.add_argument("--contact-cutoff", type=float, default=4.5,
                    help="heavy-atom A defining a residue-residue contact")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    for target in args.targets:
        jobs = discover(Path(args.boltz_out), target)
        if not jobs:
            print(f"[{target}] no predictions found")
            continue
        rows = []
        if args.workers > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futs = [pool.submit(score_one, *j, args.dist_cutoff,
                                    args.contact_cutoff) for j in jobs]
                for n, fut in enumerate(as_completed(futs), 1):
                    rows.append(fut.result())
                    if n % 500 == 0:
                        print(f"  [{target}] {n}/{len(jobs)}", flush=True)
        else:
            for j in jobs:
                rows.append(score_one(*j, args.dist_cutoff, args.contact_cutoff))
        rows.sort(key=lambda r: r["design_id"])
        out = Path(args.out_dir) / f"interface_metrics_{target}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS)
            w.writeheader()
            w.writerows(rows)
        ok = sum(1 for r in rows if r["status"] == "ok")
        print(f"[{target}] interface metrics ok={ok}/{len(rows)} -> "
              f"{os.path.relpath(out, REPO)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
