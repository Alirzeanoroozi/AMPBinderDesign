#!/usr/bin/env python3
"""Run RFdiffusion binder-backbone generation for the stage-1 targets.

RFdiffusion generates backbone complexes. The designed binder chain is poly-Gly
until a sequence-design method is run; this script writes the BoltzGen-style
metadata needed to inverse-fold those RFdiffusion backbones with
``src/1_design/inverse_folding.py``.

Example dry run:
  python src/1_design/rfdiffusion.py --rfdiffusion-dir /path/to/RFdiffusion --dry-run

Example real run:
  python src/1_design/rfdiffusion.py \
      --rfdiffusion-dir /path/to/RFdiffusion \
      --model-dir /path/to/RFdiffusion/models \
      --num-designs 1000
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

DEFAULT_TARGETS = ("NDM5", "KPC3")
DEFAULT_BINDER_LEN = "12-45"
DEFAULT_NOISE = 0.5


def load_target_defs() -> Dict[str, dict]:
    targets_json = REPO / "targets" / "targets.json"
    with targets_json.open() as fh:
        return {k: v for k, v in json.load(fh).items() if not k.startswith("_")}


def read_hotspot_resnums(target: str) -> List[int]:
    path = REPO / "targets" / target / "epitope" / "hotspots_union.csv"
    with path.open() as fh:
        rows = csv.DictReader(fh)
        return [int(row["scaffold_resnum"]) for row in rows]


def pdb_residues(path: Path, chain_id: Optional[str] = None) -> List[Tuple[str, int, str]]:
    residues: List[Tuple[str, int, str]] = []
    seen: set = set()
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip() or "_"
            if chain_id is not None and chain != chain_id:
                continue
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            key = (chain, resnum, icode)
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def contiguous_spans(resnums: Iterable[int]) -> List[Tuple[int, int]]:
    nums = sorted(set(resnums))
    if not nums:
        return []
    spans = []
    start = prev = nums[0]
    for num in nums[1:]:
        if num == prev + 1:
            prev = num
            continue
        spans.append((start, prev))
        start = prev = num
    spans.append((start, prev))
    return spans


def contig_for_target(pdb_path: Path, chain_id: str, binder_len: str) -> str:
    residues = pdb_residues(pdb_path, chain_id)
    if not residues:
        raise ValueError(f"{pdb_path}: no ATOM residues for chain {chain_id}")
    spans = contiguous_spans(resnum for _chain, resnum, icode in residues if not icode)
    if not spans:
        raise ValueError(f"{pdb_path}: no plain residue numbers for chain {chain_id}")
    target_contig = "/".join(f"{chain_id}{start}-{end}" for start, end in spans)
    return f"[{target_contig}/0 {binder_len}]"


def select_hotspots(hotspots: List[int], limit: int) -> List[int]:
    if limit <= 0 or len(hotspots) <= limit:
        return hotspots
    if limit == 1:
        return [hotspots[len(hotspots) // 2]]
    # Evenly sample the Stage-0 epitope so a small RFdiffusion hotspot set still
    # covers the target pocket instead of only the first residues in the CSV.
    last = len(hotspots) - 1
    idxs = {round(i * last / (limit - 1)) for i in range(limit)}
    return [hotspots[i] for i in sorted(idxs)]


def hotspot_arg(chain_id: str, hotspots: List[int]) -> str:
    values = ",".join(f"{chain_id}{resnum}" for resnum in hotspots)
    return f"[{values}]"


def validate_hotspots(pdb_path: Path, chain_id: str, hotspots: List[int]) -> None:
    present = {resnum for _chain, resnum, _icode in pdb_residues(pdb_path, chain_id)}
    missing = [resnum for resnum in hotspots if resnum not in present]
    if missing:
        raise ValueError(
            f"{pdb_path}: hotspot residues absent from chain {chain_id}: {missing}"
        )


def rfdiffusion_command(
    *,
    rfdiffusion_dir: Path,
    target: str,
    target_def: dict,
    out_dir: Path,
    num_designs: int,
    binder_len: str,
    hotspot_limit: int,
    noise: float,
    model_dir: Optional[Path],
    ckpt_override_path: Optional[Path],
    extra_args: List[str],
) -> Tuple[List[str], Path, List[int]]:
    chain_id = target_def.get("scaffold_chain", "A")
    input_pdb = REPO / "targets" / target / "structures" / f"{target}_target.pdb"
    if not input_pdb.exists():
        raise FileNotFoundError(input_pdb)

    all_hotspots = read_hotspot_resnums(target)
    selected_hotspots = select_hotspots(all_hotspots, hotspot_limit)
    validate_hotspots(input_pdb, chain_id, selected_hotspots)

    target_out = out_dir / target
    output_prefix = target_out / target
    command = [
        str(rfdiffusion_dir / "scripts" / "run_inference.py"),
        f"inference.output_prefix={output_prefix}",
        f"inference.input_pdb={input_pdb}",
        f"contigmap.contigs={contig_for_target(input_pdb, chain_id, binder_len)}",
        f"ppi.hotspot_res={hotspot_arg(chain_id, selected_hotspots)}",
        f"inference.num_designs={num_designs}",
        f"denoiser.noise_scale_ca={noise}",
        f"denoiser.noise_scale_frame={noise}",
    ]
    if model_dir is not None:
        command.append(f"inference.model_directory_path={model_dir}")
    if ckpt_override_path is not None:
        command.append(f"inference.ckpt_override_path={ckpt_override_path}")
    command.extend(extra_args)
    return command, target_out, selected_hotspots


def structure_design_mask(path: Path, target_chain: str) -> List[bool]:
    residues = pdb_residues(path)
    if not residues:
        raise ValueError(f"{path}: no ATOM residues found")
    chains = {chain for chain, _resnum, _icode in residues}
    if target_chain not in chains:
        raise ValueError(f"{path}: target chain {target_chain!r} not present")
    if len(chains) < 2:
        raise ValueError(f"{path}: no binder chain found")
    return [chain != target_chain for chain, _resnum, _icode in residues]


def write_inverse_folding_metadata(out_dir: Path, target: str, target_chain: str) -> int:
    import numpy as np

    wrote = 0
    for pdb_path in sorted(out_dir.glob(f"{target}_*.pdb")):
        mask = np.array(structure_design_mask(pdb_path, target_chain), dtype=np.bool_)
        np.savez(
            pdb_path.with_suffix(".npz"),
            design_mask=mask,
            inverse_fold_design_mask=mask.astype(np.float32),
        )
        wrote += 1
    return wrote


def write_manifest(rows: List[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", nargs="+", default=list(DEFAULT_TARGETS))
    ap.add_argument(
        "--rfdiffusion-dir",
        default=os.environ.get("RFDIFFUSION_REPO", ""),
        help="RFdiffusion checkout containing scripts/run_inference.py",
    )
    ap.add_argument(
        "--model-dir",
        default=os.environ.get("RFDIFFUSION_MODEL_DIR", ""),
        help="optional RFdiffusion model directory",
    )
    ap.add_argument("--out-dir", default=str(REPO / "rfdiffusion_outputs"))
    ap.add_argument("--num-designs", type=int, default=1000)
    ap.add_argument("--binder-len", default=DEFAULT_BINDER_LEN)
    ap.add_argument(
        "--hotspot-limit",
        type=int,
        default=0,
        help="0 keeps all Stage-0 hotspots; RFdiffusion often benefits from pilot runs with 3-6",
    )
    ap.add_argument("--noise", type=float, default=DEFAULT_NOISE)
    ap.add_argument("--ckpt-override-path", default="")
    ap.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="additional Hydra override passed through to RFdiffusion; repeatable",
    )
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--skip-metadata",
        action="store_true",
        help="do not write .npz design masks after RFdiffusion finishes",
    )
    args = ap.parse_args()

    if not args.rfdiffusion_dir:
        raise SystemExit("Set --rfdiffusion-dir or RFDIFFUSION_REPO")
    rfdiffusion_dir = Path(args.rfdiffusion_dir).expanduser().resolve()
    run_inference = rfdiffusion_dir / "scripts" / "run_inference.py"
    if not args.dry_run and not run_inference.exists():
        raise FileNotFoundError(run_inference)

    model_dir = Path(args.model_dir).expanduser().resolve() if args.model_dir else None
    ckpt_override_path = (
        Path(args.ckpt_override_path).expanduser().resolve()
        if args.ckpt_override_path
        else None
    )
    out_dir = Path(args.out_dir).expanduser().resolve()
    target_defs = load_target_defs()
    rows = []

    for target in args.targets:
        if target not in target_defs:
            raise ValueError(f"unknown target {target!r}; available: {sorted(target_defs)}")
        command, target_out, selected_hotspots = rfdiffusion_command(
            rfdiffusion_dir=rfdiffusion_dir,
            target=target,
            target_def=target_defs[target],
            out_dir=out_dir,
            num_designs=args.num_designs,
            binder_len=args.binder_len,
            hotspot_limit=args.hotspot_limit,
            noise=args.noise,
            model_dir=model_dir,
            ckpt_override_path=ckpt_override_path,
            extra_args=args.extra_arg,
        )
        chain_id = target_defs[target].get("scaffold_chain", "A")
        print(f"\n[{target}] {' '.join(command)}")
        rows.append(
            {
                "target": target,
                "input_pdb": os.path.relpath(
                    REPO / "targets" / target / "structures" / f"{target}_target.pdb",
                    REPO,
                ),
                "target_chain": chain_id,
                "binder_len": args.binder_len,
                "num_designs": args.num_designs,
                "hotspots": ",".join(f"{chain_id}{h}" for h in selected_hotspots),
                "output_dir": os.path.relpath(target_out, REPO),
                "command": " ".join(command),
            }
        )
        if args.dry_run:
            continue
        target_out.mkdir(parents=True, exist_ok=True)
        subprocess.run(command, cwd=str(rfdiffusion_dir), check=True)
        if not args.skip_metadata:
            wrote = write_inverse_folding_metadata(target_out, target, chain_id)
            print(f"[{target}] wrote {wrote} inverse-folding metadata files")

    manifest = REPO / "rfdiffusion_inputs" / "manifest.csv"
    write_manifest(rows, manifest)
    print(f"\nManifest -> {os.path.relpath(manifest, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
