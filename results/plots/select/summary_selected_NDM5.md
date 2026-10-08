AMPBinderDesign panel summary: NDM5 selected panel  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/selected_NDM5.csv  (8 designs)

INHIBITOR GATE  catalytic_ok=True: 8/8  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 5/8  (62.5%)

STRUCTURE GATES  pass: 5/8 (62.5%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                   8     0.284     0.421     0.539      0.87     0.879
ipSAE (max of A→B, B→A)                        8         0    0.0124    0.0779     0.431     0.579
ipSAE_min (d0-floored for short binders)       8         0         0    0.0126     0.397     0.495
Active-site coverage                           8       0.2      0.45       0.5       0.6      0.75
Interface focus on epitope                     8     0.174     0.429     0.462     0.571     0.714
# catalytic-core contacts                      8         1         1       1.5         2         4
pDockQ                                         8     0.149     0.221     0.272     0.338     0.348
LIS                                            8    0.0183    0.0662     0.105     0.358      0.41
Boltz-2 complex pLDDT                          8     0.857     0.875     0.912     0.935     0.937
AMPScanner P(AMP annotation)                   8    0.0026    0.0105    0.0343    0.0761      0.44
Macrel P(AMP annotation)                       8         0      0.01      0.03     0.089     0.089
Macrel P(hemolytic)                            8     0.079     0.089     0.104     0.178     0.198
Binder length (aa)                             8        53        58      68.5        81        96
Net charge at pH 7.4                           8     -7.02     -4.03     -2.53     -1.02     -0.03
Periplasmic-delivery proxy                     8       0.5       0.5       0.5       0.5       0.5
Aggregation proxy                              8     0.857      1.17      1.54      2.19      2.37
# synthesis liabilities                        8         0         1         1         2         2
Composite rank score                           8     0.653      0.74     0.746     0.838     0.891

TOP 8 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  NDM5_0_2                       0.891   0.87   0.431   0.40     1   53    -6.0
  NDM5_16_2                      0.838   0.88   0.579   0.45     1   77    -4.0
  NDM5_3_9                       0.771   0.54   0.124   0.60     2   79    -3.0
  NDM5_1_1                       0.749   0.57   0.060   0.55     2   60    -0.0
  NDM5_17_3                      0.743   0.41   0.012   0.60     2   96    -7.0
  rank001_NDM5_026               0.740   0.42   0.012   0.75     4   58    -1.0
  NDM5_11_5                      0.705   0.28   0.000   0.45     1   81    -2.0
  NDM5_10_1                      0.653   0.54   0.096   0.20     1   54    -2.0

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
