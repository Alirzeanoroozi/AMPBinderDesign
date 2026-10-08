AMPBinderDesign panel summary: KPC3  /scratch/esevinc22/AMPBinderDesign/AMPBinderDesign/results/ranked_KPC3.csv  (45 designs)

INHIBITOR GATE  catalytic_ok=True: 15/45  (33.3%)   [STRONG]

STRUCTURE GATE  catalytic_ok and iPTM>=0.5: 2/45  (4.4%)

STRUCTURE GATES  pass: 1/45 (2.2%)

metric                                         n       min       q25    median       q75       max
----------------------------------------------------------------------------------------------------
Boltz-2 iPTM                                  45      0.12     0.184     0.229     0.288     0.597
ipSAE (max of A→B, B→A)                       45         0         0         0         0    0.0723
ipSAE_min (d0-floored for short binders)      45         0         0         0         0    0.0244
Active-site coverage                          45         0     0.111     0.222     0.278     0.556
Interface focus on epitope                    45         0     0.095       0.2     0.263       0.5
# catalytic-core contacts                     45         0         0         0         1         1
pDockQ                                        45    0.0591      0.14     0.205     0.272     0.418
LIS                                           45         0         0         0    0.0116      0.12
Boltz-2 complex pLDDT                         45     0.877     0.894     0.902      0.91     0.934
AMPScanner P(AMP annotation)                  45    0.0006    0.0061    0.0112    0.0285     0.217
Macrel P(AMP annotation)                      45         0      0.01      0.03     0.069     0.426
Macrel P(hemolytic)                           45      0.02     0.059     0.089     0.109     0.208
Binder length (aa)                            45        57        65        76        86        99
Net charge at pH 7.4                          45       -13     -8.02     -6.02     -5.02     -2.03
Periplasmic-delivery proxy                    45     0.367       0.5       0.5       0.5       0.5
Aggregation proxy                             45     0.657      1.24       1.5      1.67       2.4
# synthesis liabilities                       45         0         1         1         2         2
Composite rank score                          45     0.298     0.474     0.556     0.656     0.837

TOP 10 BY rank_score:
  design_id                       rank   iPTM   ipSAE   eRec  #cat    L  charge
  KPC3_1_3                       0.837   0.60   0.072   0.56     1   57    -5.0
  rank003_KPC3_028               0.776   0.55   0.013   0.17     1   65    -9.0
  KPC3_6_5                       0.763   0.30   0.000   0.56     1   96   -13.0
  rank007_KPC3_020               0.747   0.41   0.000   0.44     0   65    -5.0
  KPC3_4_9                       0.747   0.37   0.011   0.33     0   84    -5.0
  KPC3_12_1                      0.698   0.28   0.000   0.50     0   88    -2.0
  KPC3_1_4                       0.660   0.25   0.000   0.28     1   57    -8.0
  rank005_KPC3_108               0.659   0.45   0.000   0.22     0   77    -7.0
  KPC3_5_5                       0.658   0.25   0.000   0.28     1   75    -6.0
  KPC3_15_6                      0.658   0.27   0.000   0.28     0   86    -9.0

WET-LAB READINESS:
  Pool is filtered for periplasmic-delivery proxy and generic peptide
  developability, not traditional AMP-likeness. Rank by Boltz-2 iPTM,
  ipSAE, and epitope coverage; require catalytic_ok for the shipped panel.
  - catalytic_ok + iPTM>=0.5 yield large enough for n=25: ship that panel.
  - ipSAE is the max of the two asymmetric directions (Dunbrack's reported
    value). ipSAE_min is also shown but is pinned near the d0 floor for
    binders <= 27 residues, so it does not discriminate.
  - Hemolysis and toxicity are exclusion risks for delivery-focused binders.
