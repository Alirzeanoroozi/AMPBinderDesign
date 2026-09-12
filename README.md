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

1. Generate target-conditioned peptide binder candidates for each target and
   hotspot set with BoltzGen and, optionally, RFDiffusion.
2. Annotate peptide candidates with classifier outputs, hemolysis/toxicity
   predictors, and delivery/developability descriptors.
3. Filter on Gram-negative outer-membrane/periplasmic-delivery proxy and
   generic peptide developability. The no-AMP branch does not hard-filter on
   traditional AMP-likeness.
4. Fold filtered target-binder complexes with Boltz-2 and compute structure
   confidence, ipSAE, and active-site overlap metrics.
5. Rank and select a diverse panel of the top 25 candidates per target.

## Sequence-Based Generation Methods

- BoltzGen
- RFDiffusion (target-PDB conditioned binder backbones; sequence design via the
  existing inverse-folding step)

## Stage 1 RFDiffusion

RFDiffusion is optional in `scripts/run_1_design.slurm`. It expects an installed
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

The RFDiffusion stage writes backbone complexes to `rfdiffusion_outputs/<TARGET>/`
and `.npz` design masks so `src/1_design/inverse_folding.py` can assign binder
sequences before `src/1_design/extract_binders.py` writes `fastas/<TARGET>.fasta`.
Do not feed raw RFDiffusion PDBs directly into the FASTA handoff unless they have
already been sequence-designed; RFDiffusion backbones use poly-Gly placeholders.

## Sequence Filter Metrics

The sequence filter does not require traditional AMP-likeness. AMPScanner,
Macrel AMP probability, HydrAMP AMP probability, and optional APEX MIC outputs
are retained as annotations when available, but hard filtering uses:

1. positive but not excessive net charge
2. amphipathicity / periplasmic-delivery proxy
3. moderate hydrophobicity and low aggregation proxy
4. low synthesis/developability liabilities
5. non-hemolytic and non-toxic classifier calls
