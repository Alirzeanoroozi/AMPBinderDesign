AMPBinderDesign panel summary: NDM5  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/ranked_NDM5.csv  (52 designs)

INHIBITOR GATE  catalytic_ok=True: 20/52  (38.5%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 8/52  (15.4%)

STRUCTURE GATES  pass: 5/52 (9.6%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  52    0.0758     0.174     0.276     0.433     0.879
ipSAE (max of A→B, B→A)                       52         0         0         0    0.0126     0.643
ipSAE_min (d0-floored for short binders)      52         0         0         0         0     0.562
Active-site coverage                          52         0         0     0.175      0.35      0.75
Interface focus on epitope                    52         0         0     0.236     0.429     0.765
# catalytic-core contacts                     52         0         0         0         1         6
pDockQ                                        52         0    0.0604     0.138     0.253     0.348
LIS                                           52         0         0    0.0086    0.0662     0.487
Boltz-2 complex pLDDT                         52     0.856     0.896     0.909     0.916     0.938
AMPScanner P(AMP annotation)                  52    0.0012    0.0038    0.0097     0.033      0.44
Macrel P(AMP annotation)                      52         0      0.01      0.02     0.059     0.356
Macrel P(hemolytic)                           52      0.05     0.079     0.099     0.139     0.267
Binder length (aa)                            52        53        77        81        90        96
Net charge at pH 7.4                          52       -12     -8.02     -5.03     -3.03      1.97
Periplasmic-delivery proxy                    52     0.348       0.5       0.5       0.5     0.664
Aggregation proxy                             52     0.186      1.11      1.48      1.79      2.37
# synthesis liabilities                       52         0         0         1         1         2
Composite rank score                          52     0.206     0.418     0.533       0.7     0.891

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  NDM5_0_2                       0.891   0.87   0.431   0.40     1   53    -6.0
  NDM5_5_9                       0.842   0.86   0.643   0.65     6   88    -8.0
  NDM5_16_2                      0.838   0.88   0.579   0.45     1   77    -4.0
  NDM5_14_2                      0.816   0.70   0.322   0.45     1   92    -9.0
  NDM5_3_9                       0.771   0.54   0.124   0.60     2   79    -3.0
  NDM5_1_1                       0.749   0.57   0.060   0.55     2   60    -0.0
  NDM5_17_3                      0.743   0.41   0.012   0.60     2   96    -7.0
  rank001_NDM5_026               0.740   0.42   0.012   0.75     4   58    -1.0
  NDM5_11_0                      0.721   0.44   0.013   0.20     1   81    -5.0
  NDM5_11_8                      0.709   0.34   0.012   0.25     1   81    -6.0

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
