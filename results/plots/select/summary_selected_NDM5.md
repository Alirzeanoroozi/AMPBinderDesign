AMPBinderDesign panel summary: NDM5 selected panel  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/selected_NDM5.csv  (25 designs)

INHIBITOR GATE  catalytic_ok=True: 25/25  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 25/25  (100.0%)

STRUCTURE GATES  pass: 25/25 (100.0%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  25     0.725     0.815     0.844     0.893     0.911
ipSAE (max of A→B, B→A)                       25     0.359     0.493     0.549     0.697       0.8
ipSAE_min (d0-floored for short binders)      25    0.0531    0.0961     0.132     0.205     0.407
Active-site coverage                          25      0.55       0.7       0.8       0.8       0.9
Interface focus on epitope                    25     0.484     0.667     0.708     0.762     0.941
# catalytic-core contacts                     25         4         5         6         6         6
pDockQ                                        25     0.131     0.168     0.242     0.302     0.436
LIS                                           25     0.418     0.499     0.556     0.614      0.69
Boltz-2 complex pLDDT                         25     0.907     0.921     0.931     0.935     0.964
AMPScanner P(AMP annotation)                  25    0.0008    0.0056    0.0101    0.0263     0.209
Macrel P(AMP annotation)                      25      0.02     0.069     0.109     0.149     0.347
Macrel P(hemolytic)                           25         0      0.01      0.03      0.05     0.099
Binder length (aa)                            25        12        14        17        25        39
Net charge at pH 7.4                          25     -5.13     -2.03     -1.02     -0.02      1.98
Periplasmic-delivery proxy                    25     0.168     0.289      0.39     0.457      0.63
Aggregation proxy                             25     -1.03       0.2       0.6      1.13      1.76
# synthesis liabilities                       25         0         0         0         1         2
Composite rank score                          25     0.889     0.894     0.905     0.922     0.933

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
