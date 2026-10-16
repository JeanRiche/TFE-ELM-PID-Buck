# Metriques des cinq methodes, banc commun

Definitions : en-tete de metriques_banc.py. Grandeurs descriptives ; le classement repose sur l'IAE de la fenetre de classement.

## Resume : IAE de classement (mV.s) ; erreur moyenne absolue en fin d'essai (mV)

| Methode | S1 | S2 | S3 | S8a | S10 |
| --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 93.77 ; 0.01 | 2.84 ; 0.00 | 8.67 ; 0.00 | 12.07 ; 0.00 | 70.48 ; 9.23 |
| PSO-PID | 87.62 ; 0.01 | 2.04 ; 0.01 | 5.46 ; 0.00 | 2.46 ; 0.00 | 3.98 ; 0.00 |
| Fuzzy-PID | 110.91 ; 0.02 | 2.43 ; 0.01 | 7.94 ; 0.03 | 4.12 ; 0.00 | 9.18 ; 0.01 |
| ELM-PID | 93.72 ; 0.00 | 2.41 ; 0.01 | 6.41 ; 0.00 | 2.96 ; 0.00 | 8.63 ; 0.01 |
| PINN-PID | 92.06 ; 0.00 | 2.44 ; 0.00 | 6.40 ; 0.01 | 2.95 ; 0.00 | 6.63 ; 0.02 |

## Classement de chaque scenario (IAE de la fenetre de classement, mV.s ; "=" : moins de 5 % ; "*" : un evenement de la fenetre non revenu)

- S1 (IAE de 0 a 30 ms) : PSO-PID 87.62 < PINN-PID 92.06 = ELM-PID 93.72 = Ziegler-Nichols 93.77 < Fuzzy-PID 110.91
- S2 (IAE de 30 ms a la fin) : PSO-PID 2.04 < ELM-PID 2.41 = Fuzzy-PID 2.43 = PINN-PID 2.44 < Ziegler-Nichols 2.84
- S3 (IAE de 30 ms a la fin) : PSO-PID 5.46 < PINN-PID 6.40 = ELM-PID 6.41 < Fuzzy-PID 7.94 < Ziegler-Nichols 8.67
- S8a (IAE de 30 ms a la fin) : PSO-PID 2.46 < PINN-PID 2.95 = ELM-PID 2.96 < Fuzzy-PID 4.12 < Ziegler-Nichols 12.07
- S10 (IAE de 100 ms a la fin) : PSO-PID 3.98 < PINN-PID 6.63 < ELM-PID 8.63 < Fuzzy-PID 9.18
  Ziegler-Nichols (reference, hors classement) 70.48 mV.s, 4 evenement(s) de la fenetre non revenu(s) ; IAE / IAE de ZN : PSO-PID 0.056, Fuzzy-PID 0.130, ELM-PID 0.122, PINN-PID 0.094

## S1 : demarrage, fenetre de classement (IAE de 0 a 30 ms), commande

| Methode | montee 10-90 % (ms) | depassement (%) | etabli +-1 V (ms) | iL crete (A) | IAE (mV.s) | ISE (V2.ms) | ITAE (mV.s2) | e max (V) | e eff (mV) | |dd| moyen | butee (%) | DCM (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.09 | 13.99 | 3.40 | 24.14 | 93.77 | 6256 | 0.0671 | 100.00 | 14440.5 | 0.0033 | 6.0 | 0.2 |
| PSO-PID | 1.09 | 7.06 | 2.64 | 22.97 | 87.62 | 6064 | 0.05479 | 100.00 | 14217.7 | 0.0080 | 5.3 | 0.2 |
| Fuzzy-PID | 1.09 | 27.34 | 3.35 | 26.52 | 110.91 | 6687 | 0.1094 | 100.00 | 14930.0 | 0.0079 | 7.5 | 0.2 |
| ELM-PID | 1.09 | 13.99 | 3.40 | 24.14 | 93.72 | 6256 | 0.06685 | 100.00 | 14440.5 | 0.0038 | 6.0 | 0.2 |
| PINN-PID | 1.09 | 13.20 | 2.70 | 24.03 | 92.06 | 6243 | 0.06219 | 100.00 | 14425.9 | 0.0036 | 5.9 | 0.2 |

## S2 : demarrage, fenetre de classement (IAE de 30 ms a la fin), commande

| Methode | montee 10-90 % (ms) | depassement (%) | etabli +-1 V (ms) | iL crete (A) | IAE (mV.s) | ISE (V2.ms) | ITAE (mV.s2) | e max (V) | e eff (mV) | |dd| moyen | butee (%) | DCM (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.09 | 13.99 | 3.40 | 24.14 | 2.84 | 1.074 | 0.1623 | 2.01 | 79.5 | 0.0032 | 0.0 | 0.0 |
| PSO-PID | 1.09 | 7.06 | 2.64 | 22.97 | 2.04 | 0.5503 | 0.1406 | 2.01 | 56.9 | 0.0082 | 0.0 | 0.0 |
| Fuzzy-PID | 1.09 | 27.34 | 3.35 | 26.52 | 2.43 | 0.8441 | 0.1498 | 2.01 | 70.5 | 0.0057 | 0.0 | 0.0 |
| ELM-PID | 1.09 | 13.99 | 3.40 | 24.14 | 2.41 | 0.7985 | 0.1508 | 2.01 | 68.5 | 0.0038 | 0.0 | 0.0 |
| PINN-PID | 1.09 | 13.20 | 2.70 | 24.03 | 2.44 | 0.8296 | 0.1517 | 2.01 | 69.9 | 0.0033 | 0.0 | 0.0 |

## S3 : demarrage, fenetre de classement (IAE de 30 ms a la fin), commande

| Methode | montee 10-90 % (ms) | depassement (%) | etabli +-1 V (ms) | iL crete (A) | IAE (mV.s) | ISE (V2.ms) | ITAE (mV.s2) | e max (V) | e eff (mV) | |dd| moyen | butee (%) | DCM (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.09 | 13.99 | 3.40 | 24.14 | 8.67 | 21.96 | 0.3494 | 5.73 | 359.4 | 0.0032 | 0.2 | 0.0 |
| PSO-PID | 1.09 | 7.06 | 2.64 | 22.97 | 5.46 | 14.95 | 0.2464 | 5.52 | 296.6 | 0.0081 | 0.4 | 0.0 |
| Fuzzy-PID | 1.09 | 27.34 | 3.35 | 26.52 | 7.94 | 20.24 | 0.325 | 5.66 | 345.1 | 0.0057 | 0.2 | 0.0 |
| ELM-PID | 1.09 | 13.99 | 3.40 | 24.14 | 6.41 | 17.06 | 0.2774 | 5.66 | 316.7 | 0.0039 | 0.3 | 0.0 |
| PINN-PID | 1.09 | 13.20 | 2.70 | 24.03 | 6.40 | 17.1 | 0.2758 | 5.66 | 317.1 | 0.0037 | 0.3 | 0.0 |

## S8a : demarrage, fenetre de classement (IAE de 30 ms a la fin), commande

| Methode | montee 10-90 % (ms) | depassement (%) | etabli +-1 V (ms) | iL crete (A) | IAE (mV.s) | ISE (V2.ms) | ITAE (mV.s2) | e max (V) | e eff (mV) | |dd| moyen | butee (%) | DCM (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 0.54 | 52.23 | 7.44 | 13.15 | 12.07 | 5.123 | 0.6202 | 1.43 | 173.6 | 0.0031 | 0.0 | 0.0 |
| PSO-PID | 0.54 | 43.49 | 3.55 | 12.60 | 2.46 | 0.3174 | 0.1784 | 0.60 | 43.2 | 0.0081 | 0.0 | 0.0 |
| Fuzzy-PID | 0.54 | 58.70 | 3.15 | 13.20 | 4.12 | 1.97 | 0.2638 | 1.40 | 107.6 | 0.0057 | 0.0 | 0.0 |
| ELM-PID | 0.54 | 52.23 | 5.36 | 13.15 | 2.96 | 0.6064 | 0.1902 | 0.97 | 59.7 | 0.0048 | 0.0 | 0.0 |
| PINN-PID | 0.54 | 52.14 | 4.28 | 13.14 | 2.95 | 0.726 | 0.2 | 0.94 | 65.3 | 0.0058 | 0.0 | 0.0 |

## S10 : demarrage, fenetre de classement (IAE de 100 ms a la fin), commande

| Methode | montee 10-90 % (ms) | depassement (%) | etabli +-1 V (ms) | iL crete (A) | IAE (mV.s) | ISE (V2.ms) | ITAE (mV.s2) | e max (V) | e eff (mV) | |dd| moyen | butee (%) | DCM (%) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.09 | 13.99 | 3.40 | 24.14 | 70.48 | 92.86 | 3.103 | 3.22 | 963.6 | 0.0029 | 0.0 | 0.0 |
| PSO-PID | 1.09 | 7.06 | 2.64 | 22.97 | 3.98 | 2.976 | 0.1572 | 1.87 | 172.5 | 0.0062 | 0.4 | 0.0 |
| Fuzzy-PID | 1.09 | 27.34 | 3.35 | 26.52 | 9.18 | 9.866 | 0.3596 | 2.18 | 314.1 | 0.0045 | 0.0 | 0.0 |
| ELM-PID | 1.09 | 13.99 | 3.40 | 24.14 | 8.63 | 6.468 | 0.3088 | 1.91 | 254.3 | 0.0041 | 0.6 | 0.0 |
| PINN-PID | 1.09 | 13.20 | 2.70 | 24.03 | 6.63 | 5.958 | 0.2605 | 2.05 | 244.1 | 0.0034 | 0.1 | 0.0 |

## S1 : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)

| Fenetre | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 10 dernieres ms | +0.01 / 0.01 / 27.0 | +0.01 / 0.01 / 27.0 | -0.02 / 0.02 / 27.5 | -0.00 / 0.00 / 27.0 | -0.00 / 0.00 / 26.9 |

## S2 : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)

| Fenetre | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 5 ms avant 50 ms | +0.00 / 0.00 / 27.0 | -0.01 / 0.01 / 27.0 | +0.01 / 0.01 / 27.6 | -0.00 / 0.00 / 27.0 | +0.02 / 0.02 / 26.9 |
| 5 ms avant 70 ms | -0.37 / 0.37 / 185.4 | -0.32 / 0.32 / 183.5 | -0.52 / 0.52 / 185.3 | -0.29 / 0.29 / 182.2 | -0.38 / 0.38 / 183.5 |
| 10 dernieres ms | -0.00 / 0.00 / 27.1 | -0.01 / 0.01 / 27.0 | -0.01 / 0.01 / 27.5 | +0.01 / 0.01 / 27.0 | +0.00 / 0.00 / 27.0 |

## S2 : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V

| Evenement | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 50 ms | 1.14 / 2.01 / 0.60 ms | 0.59 / 2.01 / 0.13 ms | 0.94 / 2.01 / 0.20 ms | 0.82 / 2.01 / 0.16 ms | 0.84 / 2.01 / 0.17 ms |
| 70 ms | 1.53 / 1.11 / 0.03 ms | 1.28 / 1.11 / 0.03 ms | 1.32 / 1.11 / 0.03 ms | 1.41 / 1.11 / 0.03 ms | 1.42 / 1.11 / 0.03 ms |

## S3 : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)

| Fenetre | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 5 ms avant 50 ms | +0.00 / 0.00 / 27.0 | -0.01 / 0.01 / 27.0 | +0.01 / 0.01 / 27.6 | -0.00 / 0.00 / 27.0 | +0.02 / 0.02 / 26.9 |
| 5 ms avant 70 ms | -0.00 / 0.00 / 27.0 | +0.00 / 0.00 / 26.9 | +0.01 / 0.01 / 27.1 | +0.01 / 0.01 / 26.9 | -0.00 / 0.00 / 26.9 |
| 10 dernieres ms | +0.00 / 0.00 / 27.0 | -0.00 / 0.00 / 27.0 | +0.03 / 0.03 / 27.5 | +0.00 / 0.00 / 27.1 | -0.01 / 0.01 / 27.0 |

## S3 : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V

| Evenement | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 50 ms | 3.64 / 5.52 / 1.05 ms | 2.20 / 5.52 / 0.46 ms | 3.36 / 5.52 / 1.31 ms | 2.57 / 5.52 / 0.87 ms | 2.64 / 5.52 / 0.91 ms |
| 70 ms | 4.85 / 5.73 / 1.49 ms | 3.09 / 5.43 / 0.45 ms | 4.40 / 5.66 / 1.34 ms | 3.66 / 5.66 / 0.87 ms | 3.58 / 5.66 / 0.85 ms |

## S8a : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)

| Fenetre | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 5 ms avant 50 ms | +0.01 / 0.01 / 27.3 | -0.01 / 0.01 / 27.3 | +0.02 / 0.02 / 27.3 | +0.00 / 0.00 / 27.4 | -0.00 / 0.00 / 27.2 |
| 5 ms avant 70 ms | +3.09 / 3.09 / 252.4 | -0.02 / 0.02 / 21.6 | +0.03 / 0.03 / 21.5 | +0.02 / 0.02 / 22.7 | +0.00 / 0.00 / 21.6 |
| 5 ms avant 100 ms | +0.01 / 0.01 / 27.5 | -0.01 / 0.01 / 27.3 | +0.02 / 0.02 / 27.5 | +0.01 / 0.01 / 28.1 | -0.01 / 0.01 / 27.3 |
| 5 ms avant 120 ms | +0.02 / 0.02 / 33.7 | -0.00 / 0.00 / 33.6 | -0.02 / 0.02 / 33.5 | +0.08 / 0.08 / 35.4 | -0.00 / 0.00 / 33.6 |
| 10 dernieres ms | -0.00 / 0.00 / 27.3 | -0.00 / 0.00 / 27.3 | -0.00 / 0.00 / 27.5 | -0.00 / 0.00 / 27.8 | +0.00 / 0.00 / 27.3 |

## S8a : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V

| Evenement | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 50 ms | 5.22 / 1.35 / 0.57 ms | 0.42 / 0.56 / dans la bande | 0.91 / 1.14 / 0.61 ms | 0.83 / 0.97 / dans la bande | 0.64 / 0.94 / dans la bande |
| 70 ms | 3.01 / 1.43 / 0.53 ms | 0.55 / 0.60 / dans la bande | 1.06 / 1.40 / 0.62 ms | 0.67 / 0.77 / dans la bande | 0.68 / 0.91 / dans la bande |
| 100 ms | 1.09 / 0.98 / dans la bande | 0.41 / 0.44 / dans la bande | 0.68 / 0.98 / dans la bande | 0.39 / 0.49 / dans la bande | 0.48 / 0.60 / dans la bande |
| 120 ms | 2.59 / 0.93 / dans la bande | 0.90 / 0.40 / dans la bande | 1.29 / 0.76 / dans la bande | 0.90 / 0.42 / dans la bande | 0.98 / 0.52 / dans la bande |

## S10 : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)

| Fenetre | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 5 ms avant 50 ms | +0.00 / 0.00 / 27.0 | -0.01 / 0.01 / 27.0 | +0.01 / 0.01 / 27.6 | -0.00 / 0.00 / 27.0 | +0.02 / 0.02 / 26.9 |
| 5 ms avant 100 ms | +0.51 / 0.51 / 42.7 | +0.02 / 0.02 / 21.8 | +0.02 / 0.02 / 21.5 | +0.02 / 0.02 / 23.2 | +0.03 / 0.03 / 21.8 |
| 5 ms avant 115 ms | -3.39 / 3.39 / 788.5 | +0.00 / 0.00 / 21.1 | +0.01 / 0.01 / 21.2 | +0.00 / 0.00 / 21.4 | +0.00 / 0.00 / 21.1 |
| 5 ms avant 130 ms | +35.65 / 35.65 / 1751.8 | +0.02 / 0.02 / 21.6 | +0.00 / 0.00 / 21.5 | +0.08 / 0.08 / 22.6 | +0.00 / 0.00 / 21.4 |
| 5 ms avant 145 ms | +63.33 / 63.33 / 2522.7 | +0.01 / 0.01 / 21.3 | +0.00 / 0.00 / 21.4 | -0.02 / 0.02 / 22.1 | -0.01 / 0.01 / 21.3 |
| 5 ms avant 160 ms | -13.63 / 13.63 / 733.2 | +0.01 / 0.01 / 21.6 | +0.01 / 0.01 / 21.4 | +0.03 / 0.03 / 23.5 | -0.00 / 0.00 / 21.8 |
| 5 ms avant 175 ms | -3.82 / 3.82 / 820.2 | -0.01 / 0.01 / 21.1 | +0.02 / 0.02 / 21.3 | -0.00 / 0.00 / 21.4 | -0.12 / 0.12 / 21.5 |
| 10 dernieres ms | -9.23 / 9.23 / 968.0 | -0.00 / 0.00 / 21.7 | +0.01 / 0.01 / 21.5 | +0.01 / 0.01 / 23.5 | +0.02 / 0.02 / 21.5 |

## S10 : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V

| Evenement | Ziegler-Nichols | PSO-PID | Fuzzy-PID | ELM-PID | PINN-PID |
| --- | --- | --- | --- | --- | --- |
| 50 ms | 149.73 / 110.59 / 14.90 ms | 124.86 / 110.58 / 3.30 ms | 122.74 / 110.59 / 2.86 ms | 129.89 / 110.58 / 4.82 ms | 128.05 / 110.59 / 3.86 ms |
| 100 ms | 9.51 / 3.08 / pas revenu | 0.78 / 1.87 / 0.30 ms | 1.66 / 2.11 / 1.00 ms | 1.71 / 1.87 / 0.90 ms | 1.17 / 1.92 / 0.70 ms |
| 115 ms | 14.07 / 3.20 / pas revenu | 0.69 / 1.44 / 0.22 ms | 1.72 / 2.18 / 1.00 ms | 2.03 / 1.91 / 1.22 ms | 1.28 / 1.96 / 0.73 ms |
| 130 ms | 15.49 / 2.62 / pas revenu | 0.50 / 0.98 / dans la bande | 1.18 / 1.46 / 0.85 ms | 1.22 / 1.21 / 0.53 ms | 0.90 / 1.36 / 0.27 ms |
| 145 ms | 5.62 / 1.35 / 1.83 ms | 0.47 / 0.96 / dans la bande | 1.16 / 1.38 / 0.88 ms | 0.94 / 1.06 / 0.16 ms | 0.81 / 1.30 / 0.26 ms |
| 160 ms | 9.97 / 3.22 / pas revenu | 0.78 / 1.87 / 0.30 ms | 1.66 / 2.11 / 1.00 ms | 1.32 / 1.87 / 0.55 ms | 1.14 / 1.97 / 0.69 ms |
| 175 ms | 15.82 / 3.19 / 9.39 ms | 0.75 / 1.44 / 0.22 ms | 1.79 / 2.18 / 1.00 ms | 1.42 / 1.62 / 0.55 ms | 1.33 / 2.05 / 0.71 ms |

## Cout de calcul du regulateur seul, banc Python (us par appel de reg.pas ; cinq scenarios reunis)

Ordre de grandeur RELATIF entre methodes : Python interprete, un seul processus, pas un materiel embarque. Le maximum inclut les aleas du systeme.

| Methode | appel moyen | appel max | ordinaire moyen | fin de fenetre moyen (nombre) | fin de fenetre max | avec adaptation moyen (nombre) | avec adaptation max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.39 | 3763.2 | - | - | - | - | - |
| PSO-PID | 1.43 | 27109.5 | - | - | - | - | - |
| Fuzzy-PID | 12.74 | 46031.4 | - | - | - | - | - |
| ELM-PID | 3.59 | 3876.1 | 3.22 | 43.12 (2000) | 3876.1 | 840.38 (20) | 3876.1 |
| PINN-PID | 38.93 | 26825.0 | 31.28 | 872.71 (2000) | 26825.0 | 16618.37 (100) | 26825.0 |

## Cout moyen d'un appel de reg.pas par scenario (us)

| Methode | S1 | S2 | S3 | S8a | S10 |
| --- | --- | --- | --- | --- | --- |
| Ziegler-Nichols | 1.21 | 1.49 | 1.36 | 1.56 | 1.32 |
| PSO-PID | 1.19 | 1.25 | 2.02 | 1.34 | 1.34 |
| Fuzzy-PID | 15.05 | 12.67 | 12.28 | 11.88 | 11.81 |
| ELM-PID | 3.51 | 3.51 | 3.53 | 3.54 | 3.84 |
| PINN-PID | 33.83 | 36.90 | 37.60 | 39.08 | 47.23 |
