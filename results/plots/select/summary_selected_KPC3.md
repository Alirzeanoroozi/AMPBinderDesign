AMPBinderDesign panel summary: KPC3 selected panel  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/selected_KPC3.csv  (1 designs)

INHIBITOR GATE  catalytic_ok=True: 1/1  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 1/1  (100.0%)

STRUCTURE GATES  pass: 1/1 (100.0%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------

TOP 1 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  KPC3_1_3                       0.837   0.60   0.072   0.56     1   57    -5.0

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
