AMPBinderDesign panel summary: NDM5  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/ranked_NDM5.csv  (7854 designs)

INHIBITOR GATE  catalytic_ok=True: 1718/7854  (21.9%)   [MARGINAL]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 1118/7854  (14.2%)

STRUCTURE GATES  pass: 1118/7854 (14.2%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                3187    0.0526     0.363      0.52     0.667     0.951
ipSAE (max of A→B, B→A)                     3187         0         0    0.0148     0.187     0.866
ipSAE_min (d0-floored for short binders)    3187         0         0    0.0122    0.0294      0.65
Active-site coverage                        3187         0         0      0.25       0.6       0.9
Interface focus on epitope                  3187         0         0     0.333     0.565         1
# catalytic-core contacts                   3187         0         0         1         3         6
pDockQ                                      3187         0    0.0945     0.138       0.2     0.582
LIS                                         3187         0    0.0273     0.107     0.257     0.723
Boltz-2 complex pLDDT                       3187     0.816     0.893     0.907      0.92     0.966
AMPScanner P(AMP annotation)                7854    0.0002    0.0041    0.0121    0.0405         1
Macrel P(AMP annotation)                    7854         0     0.089     0.149     0.228     0.713
Macrel P(hemolytic)                         7854         0      0.05     0.089     0.139     0.495
Binder length (aa)                          7854        12        20        29        38        45
Net charge at pH 7.4                        7854       -12     -4.02     -2.02     -0.99      5.97
Periplasmic-delivery proxy                  7854     0.021     0.373     0.481       0.5     0.998
Aggregation proxy                           7854     -2.41     0.443     0.957      1.47       2.5
# synthesis liabilities                     7854         0         0         0         1         2
Composite rank score                        7854     0.102     0.166     0.197     0.479     0.933

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  NDM5_0833                      0.933   0.91   0.733   0.80   6.0   22    -1.1
  NDM5_4051                      0.932   0.81   0.506   0.90   6.0   19    +0.9
  NDM5_4819                      0.929   0.83   0.493   0.75   4.0   13    -1.0
  NDM5_10_1                      0.927   0.81   0.390   0.80   6.0   13    -0.0
  NDM5_1778                      0.924   0.89   0.716   0.70   6.0   26    -1.0
  NDM5_423_7                     0.924   0.86   0.493   0.70   5.0   12    -1.0
  rank0070_NDM5_0374             0.922   0.88   0.721   0.80   5.0   25    -5.1
  NDM5_663_6                     0.921   0.91   0.672   0.60   6.0   13    -0.0
  NDM5_4423                      0.915   0.77   0.456   0.85   6.0   23    -3.1
  NDM5_1045                      0.910   0.83   0.537   0.80   6.0   25    -4.1

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
