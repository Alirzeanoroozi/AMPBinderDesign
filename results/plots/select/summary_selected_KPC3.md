AMPBinderDesign panel summary: KPC3 selected panel  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/selected_KPC3.csv  (25 designs)

INHIBITOR GATE  catalytic_ok=True: 25/25  (100.0%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 25/25  (100.0%)

STRUCTURE GATES  pass: 25/25 (100.0%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  25     0.764     0.832     0.852     0.881     0.904
ipSAE (max of A→B, B→A)                       25     0.381     0.524     0.602     0.639     0.739
ipSAE_min (d0-floored for short binders)      25    0.0429    0.0761     0.112      0.16     0.195
Active-site coverage                          25       0.5     0.778     0.833     0.833     0.944
Interface focus on epitope                    25     0.333     0.516     0.552       0.6      0.75
# catalytic-core contacts                     25         2         3         3         3         4
pDockQ                                        25     0.136     0.235     0.278     0.333     0.439
LIS                                           25     0.383     0.495     0.519     0.586     0.661
Boltz-2 complex pLDDT                         25     0.884     0.891     0.912     0.918      0.94
AMPScanner P(AMP annotation)                  25    0.0005    0.0022    0.0088    0.0197     0.123
Macrel P(AMP annotation)                      25      0.04     0.069     0.109     0.158     0.347
Macrel P(hemolytic)                           25         0         0      0.01      0.04     0.099
Binder length (aa)                            25        12        17        22        24        36
Net charge at pH 7.4                          25     -5.02     -2.02     -1.03     -0.03      1.97
Periplasmic-delivery proxy                    25     0.049     0.282     0.353     0.411       0.5
Aggregation proxy                             25    -0.714     0.143     0.614     0.971      2.09
# synthesis liabilities                       25         0         0         0         1         2
Composite rank score                          25     0.896     0.902     0.911     0.924     0.957

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  KPC3_3781                      0.957   0.90   0.700   0.94   3.0   18    -0.0
  KPC3_2683                      0.955   0.90   0.602   0.94   4.0   17    -5.0
  KPC3_0391                      0.937   0.89   0.739   0.78   3.0   26    -2.0
  KPC3_0472                      0.934   0.84   0.627   0.89   3.0   25    -1.0
  KPC3_236_3                     0.932   0.82   0.595   0.83   3.0   22    +1.0
  KPC3_0834                      0.928   0.90   0.703   0.83   3.0   23    -4.0
  KPC3_1228                      0.924   0.84   0.604   0.89   4.0   22    -0.0
  KPC3_1641                      0.923   0.89   0.639   0.72   2.0   13    -0.0
  KPC3_685_0                     0.921   0.90   0.734   0.83   3.0   15    +2.0
  KPC3_1052                      0.919   0.84   0.548   0.83   3.0   19    -1.0

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
