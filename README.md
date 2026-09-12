# AMPBinderDesign — peptide inhibitors for NDM-5 and KPC-3

AMPBinderDesign is a sequence-driven pipeline for designing short peptide binders against antimicrobial-resistance proteins.
The first targets are the carbapenemases NDM-5 and KPC-3, which contribute to carbapenem resistance in Klebsiella pneumoniae.

>NDM-5
MELPNIMHPVAKLSTALAAALMLSGCMPGEIRPTIGQQMETGDQRFGDLVFRQLAPNVWQHTSYLDMPGFGAVASNGLIVRDGGRVLLVDTAWTDDQTAQILNWIKQEINLPVALAVVTHAHQDKMGGMDALHAAGIATYANALSNQLAPQEGLVAAQHSLTFAANGWVEPATAPNFGPLKVFYPGPGHTSDNITVGIDGTDIAFGGCLIKDSKAKSLGNLGDADTEHYAASARAFGAAFPKASMIVMSHSAPDSRAAITHTARMADKLR

Close to 4EYL

>KPC-3
MSLYRRLVLLSCLSWPLAGFSATALTNLVAEPFAKLEQDFGGSIGVYAMDTGSGATVSYRAEERFPLCSSFKGFLAAAVLARSQQQAGLLDTPIRYGKNALVPWSPISEKYLTTGMTVAELSAAAVQYSDNAAANLLLKELGGPAGLTAFMRSIGDTTFRLDRWELELNSAIPGDARDTSSPRAVTESLQKLTLGSALAAPQRQQFVDWLKGNTTGNHRIRAAVPADWAVGDKTGTCGVYGTANDYAVVWPTGRAPIVLAVYTRAPNKDDKYSEAVIAAAARLALEGLGVNGQ

Close to 3DW0

## Pipeline

1. run Boltzgen on each sequence as target and generate binders with specific hotspots.
2. run AMPScanner, Macrel, ToxinPred, Apex and Hydramp.
3. Filter based on Boltz2 Strucutre prediction.
4. filter 25 best for each.

## Sequence-Based Generation Methods

- Boltzgen
- RFdiffusion (target-PDB conditioned binder backbones; sequence design via the
  existing inverse-folding step)

## Stage 1 RFDiffusion

RFdiffusion is optional in `scripts/run_1_design.slurm`. It expects an installed
checkout from `https://github.com/RosettaCommons/RFdiffusion.git` with model
weights available in `models/`.

Dry-run the generated commands:

```bash
python src/1_design/rfdiffusion.py \
  --rfdiffusion-dir /path/to/RFdiffusion \
  --model-dir /path/to/RFdiffusion/models \
  --dry-run
```

Run it through the stage-1 SLURM job:

```bash
RUN_RFDIFFUSION=1 \
RFDIFFUSION_REPO=/path/to/RFdiffusion \
RFDIFFUSION_MODEL_DIR=/path/to/RFdiffusion/models \
sbatch scripts/run_1_design.slurm
```

The RFdiffusion stage writes backbone complexes to `rfdiffusion_outputs/<TARGET>/`
and `.npz` design masks so `src/1_design/inverse_folding.py` can assign binder
sequences before `src/1_design/extract_binders.py` writes `fastas/<TARGET>.fasta`.
Do not feed raw RFdiffusion PDBs directly into the FASTA handoff unless they have
already been sequence-designed; RFdiffusion backbones use poly-Gly placeholders.

## Filtering Metrics other that Boltz

1. APEX
2. HydrAMP
3. Macrel
4. AMPScanner
5. ToxinPred
