AMPBinderDesign panel summary: KPC3  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/ranked_KPC3.csv  (9284 designs)

INHIBITOR GATE  catalytic_ok=True: 2018/9284  (21.7%)   [MARGINAL]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 1274/9284  (13.7%)

STRUCTURE GATES  pass: 1274/9284 (13.7%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                3188    0.0905     0.344     0.514     0.685     0.973
ipSAE (max of A→B, B→A)                     3188         0         0    0.0132     0.222     0.893
ipSAE_min (d0-floored for short binders)    3188         0         0    0.0115    0.0247     0.602
Active-site coverage                        3188         0     0.056     0.167     0.389         1
Interface focus on epitope                  3188         0     0.077       0.2     0.375         1
# catalytic-core contacts                   3188         0         0         1         1         4
pDockQ                                      3188         0    0.0978     0.137     0.196     0.629
LIS                                         3188         0    0.0128     0.081     0.245     0.725
Boltz-2 complex pLDDT                       3188      0.84     0.888     0.902     0.914     0.962
AMPScanner P(AMP annotation)                9284    0.0001    0.0036    0.0106    0.0373         1
Macrel P(AMP annotation)                    9284         0     0.099     0.158     0.228     0.752
Macrel P(hemolytic)                         9284         0      0.04     0.079     0.129     0.495
Binder length (aa)                          9284        12        19        27        37        45
Net charge at pH 7.4                        9284       -12     -4.02     -2.02     -1.02      5.97
Periplasmic-delivery proxy                  9284     0.014     0.348     0.451       0.5     0.973
Aggregation proxy                           9284     -2.21       0.4     0.957      1.49       2.5
# synthesis liabilities                     9284         0         0         1         1         2
Composite rank score                        9284    0.0991     0.162      0.19     0.422     0.957

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
