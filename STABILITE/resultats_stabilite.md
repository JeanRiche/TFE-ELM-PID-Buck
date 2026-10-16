# Stabilité linéaire de la boucle à gains figés

Produit par `stabilite.py`.

## 1. Modèle

- Convertisseur : modèle moyen sans pertes, conduction continue, état [iL ; vC], entrée d, sortie vC.
- Discrétisation : bloqueur d'ordre zéro à Tc = 1/220 000 s. Justification (banc `banc_commun.py`, `simuler`) : d(k) est calculé à l'instant k·Tc à partir de v(k·Tc) et tenu pendant les 120 pas qui suivent.
- Retard de calcul : aucun dans le modèle de référence (le banc applique d(k) à l'instant de la mesure). Variante de sensibilité : un pas Tc (z⁻¹).
- MLI : gain unité (d comparé à une dent de scie normalisée). Bloc PID commun (éq. 3.4), N = 64 122,9 rad/s, saturation et anti-emballement ignorés.
- Points : 9 coins R = 4, 5, 6 Ω × Vin = 160, 200, 240 V ; complément hors critère : 25 et 98 Ω.

## 2. Convertisseur seul

A = [[0, −1/L], [1/C, −1/(RC)]], B = [Vin/L ; 0], C = [0, 1], D = 0. À 5 Ω, 200 V : A = [[0 ; −100], [21 276,6 ; −4 255,3]], B = [20 000 ; 0].

Pôles à 5 Ω : **−578,70** et **−3 676,62 s⁻¹**. ωn = 1 458,6 rad/s, ζ = 1,459.

| R (Ω) | ζ | Pôles (s⁻¹) | Pôles discrets (module) | écart exp(pTc) |
|---|---|---|---|---|
| 4 | 1,823 | −4 883,46 ; −435,69 | 0,97805 ; 0,99802 | 0,0e+00 |
| 5 | 1,459 | −3 676,62 ; −578,70 | 0,98343 ; 0,99737 | 1,1e-16 |
| 6 | 1,216 | −2 781,04 ; −765,06 | 0,98744 ; 0,99653 | 1,1e-16 |
| 25 | 0,292 | −425,53 ± 1 395,20 j | 0,99807 ; 0,99807 | 8,7e-18 |
| 98 | 0,074 | −108,55 ± 1 454,60 j | 0,99951 ; 0,99951 | 2,6e-18 |

Les pôles ne dépendent pas de Vin (gain seul).

## 3. Boucle fermée avec les gains de Ziegler-Nichols

| R (Ω) | Vin (V) | max |z| | Marge de phase | Coupure (Hz) | Marge de gain (dB) | à (Hz) | MG par les pôles (dB) | min |1+L| | MP avec retard Tc |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 160 | 0,99366 | 34,16° | 775 | 42,87 | 26 286 | 42,87 | 0,585 | 32,89° |
| 4 | 200 | 0,99219 | 37,00° | 886 | 40,93 | 26 286 | 40,93 | 0,632 | 35,55° |
| 4 | 240 | 0,99068 | 39,53° | 993 | 39,35 | 26 286 | 39,35 | 0,672 | 37,91° |
| 5 | 160 | 0,99437 | 29,15° | 828 | 42,73 | 26 055 | 42,73 | 0,503 | 27,80° |
| 5 | 200 | 0,99300 | 31,99° | 940 | 40,79 | 26 055 | 40,79 | 0,550 | 30,45° |
| 5 | 240 | 0,99156 | 34,57° | 1 047 | 39,20 | 26 055 | 39,20 | 0,593 | 32,86° |
| 6 | 160 | 0,99514 | 24,93° | 860 | 42,63 | 25 899 | 42,63 | 0,432 | 23,52° |
| 6 | 200 | 0,99383 | 27,94° | 971 | 40,69 | 25 899 | 40,69 | 0,483 | 26,35° |
| 6 | 240 | 0,99245 | 30,68° | 1 078 | 39,11 | 25 899 | 39,11 | 0,529 | 28,91° |

La marge de gain est recoupée de façon indépendante : facteur de gain qui amène un pôle sur le cercle unité (colonne « par les pôles »). Avec le retard d'un pas, la marge de gain ZN tombe entre 29,98 et 33,78 dB et le module max reste sous 0,99534.

Charge légère (hors critère de l'ensemble, pour mémoire) :

| R (Ω) | Vin (V) | max |z| | Marge de phase | Coupure (Hz) | Marge de gain à la hausse (dB) | Marge de gain à la baisse (dB) |
|---|---|---|---|---|---|---|
| 25 | 160 | 0,99940 | 3,14° | 938 | 42,24 | 1,32 |
| 25 | 200 | 0,99825 | 8,06° | 1 046 | 40,30 | 3,26 |
| 25 | 240 | 0,99701 | 12,32° | 1 150 | 38,72 | 4,85 |
| 98 | 160 | 1,00059 | −3,09° | 943 | instable | instable |
| 98 | 200 | 0,99946 | 2,51° | 1 051 | 40,21 | 0,86 |
| 98 | 240 | 0,99824 | 7,30° | 1 155 | 38,63 | 2,44 |

À charge légère, la résonance du filtre (ζ = 0,29 et 0,074, f0 = 232 Hz) fait passer la phase de L par −180° à une fréquence où |L| > 1 : la boucle est conditionnellement stable, et une baisse du gain de quelques dB la déstabilise (marge à la baisse vérifiée par les pôles).

## 4. ELM-PID à gains figés sur l'ensemble admissible

Ensemble tel que la loi le voit : g interpolé (trilinéaire en log2) dans `ensemble_gains_elm.mat`, g ≤ 0, x ∈ [1/4 ; 4]³. Pire cas sur les 9 coins pour chaque x.

| Jeu de x | Nombre | max |z| | Pire marge de phase | Plus haute coupure (Hz) | Pire marge de gain (dB) | min |1+L| |
|---|---|---|---|---|---|---|
| Points admissibles de la table (pas 2^(1/8)) | 12663 | 0,99913 | 24,930° | 2 199,9 | 28,59 | 0,333 |
| Grille 2^(1/16), admis | 99033 | 0,99913 | 24,843° | 2 200,0 | 28,59 | 0,333 |
| Grille 2^(1/32), admis (pôles seuls) | 782796 | 0,99915 | - | - | - | - |
| Gains atteints par l'ELM-PID (S1 à S9, S10) | 45 | 0,99514 | 24,930° | 2 134,9 | 32,49 | - |

- Aucun x admis n'est instable aux 9 coins : 0 sur 782796 (grille 2^(1/32)). Pire module 0,99915 en x = (0,250 ; 0,723 ; 3,668), coin 4 Ω / 160 V : un pôle lent (petit xP), pas une oscillation.
- Pire marge de phase sur la table : 24,930° au départ (1 ; 1 ; 1), coin 6 Ω / 160 V, égale à M0 = 24,930° du JSON.
- Entre les points de la table, 23 points de la grille 2^(1/16) admis par l'interpolation ont une marge sous M0 ; la pire vaut 24,843° (écart −0,087°), en x = (0,771 ; 1,682 ; 2,484). Tous sont à la frontière (g interpolé entre -7.0e-03 et -9.4e-06). Pire marge à l'intérieur strict : 25,25°.
- Coupure : au plus 2 200,03 Hz (borne 2 200 Hz ; dépassement de 0,03 Hz, effet de l'interpolation).
- Marge de gain à la hausse : au moins 28,59 dB (x = (1,682 ; 4,000 ; 4,000)), aucun point sous 6 dB ; 77,2 % de l'ensemble est sous la plus petite marge de ZN (39,11 dB). 1743 réglages (sur 99033) sont conditionnellement stables à au moins un coin (marge à la baisse au plus 4,97 dB, x = (0,917 ; 3,084 ; 4,000)).
- Gains atteints : x entre (1,000 ; 0,864 ; 1,000) et (3,467 ; 2,417 ; 1,499), tous admissibles (g max 1.2e-13), marge de phase au moins 24,93°, marge de gain au moins 32,49 dB.

Variante avec un retard d'un pas Tc (sensibilité) : aucun x admis instable (0), pire module 0,99915 ; pire marge de phase 21,46° (x = (3,668 ; 0,545 ; 0,917)), perte de 0,41 à 3,60° selon x ; pire marge de gain 19,58 dB.

Charge légère (hors critère ; l'ensemble ne promet rien ici), grille 2^(1/16) :

| R (Ω) | Vin (V) | Part de l'ensemble stable | Gains atteints instables | max |z| ZN |
|---|---|---|---|---|
| 25 | 160 | 98,2 % | aucun | 0,99940 |
| 25 | 200 | 99,3 % | aucun | 0,99825 |
| 25 | 240 | 99,8 % | aucun | 0,99701 |
| 98 | 160 | 90,7 % | (1,000 ; 1,000 ; 1,000), (1,143 ; 1,104 ; 1,017), (1,256 ; 1,190 ; 1,040) | 1,00059 |
| 98 | 200 | 95,0 % | aucun | 0,99946 |
| 98 | 240 | 97,2 % | aucun | 0,99824 |

À 98 Ω sous 160 V, ZN lui-même est instable ; les deux autres réglages instables, voisins de ZN, ont été atteints sur S8b (98 Ω). Le calcul à gains figés ne dit pas si l'essai les a tenus pendant une phase à 160 V.

## 5. Limites : ce qui est prouvé et ce qui ne l'est pas

Prouvé, sous les hypothèses du modèle (moyen, linéaire, conduction continue, sans retard de calcul) :

- la boucle ZN est stable aux 9 coins nominaux, avec 24,9 à 39,5° de marge de phase ;
- tout réglage figé de l'ensemble admissible (table interpolée) est stable aux 9 coins, avec au moins 24,84° de marge de phase et 28,59 dB de marge de gain ; il le reste avec un retard de calcul d'un pas.

Non prouvé :

- la stabilité du système à gains variables : les gains changent toutes les 0,5 ms ; un système linéaire dont chaque gel est stable peut être instable si ses paramètres varient assez vite. Le calcul ne borne pas la vitesse d'adaptation (voir le chapitre 4 du mémoire) ;
- la saturation du rapport cyclique, l'anti-emballement, la conduction discontinue, l'échantillonnage de l'ondulation de découpage (10 échantillons par période) et la modulation multiple dans une période ;
- la stabilité à charge légère : à 98 Ω sous 160 V, ZN et 9,3 % de l'ensemble sont instables ; l'ensemble admissible ne contient aucune garantie hors de 4 à 6 Ω ;
- les points entre les coins (R et Vin intermédiaires) et les tolérances de L et C ;
- toute valeur au-delà d'environ fs/2 (11 kHz) : la marge de gain, mesurée vers 26 kHz, en fait partie.

## 6. Figures (PDF vectoriels, `figures/`)

- `bode_boucle_ouverte_zn_nominal.pdf` : Bode de L(z) ZN à 5 Ω, 200 V, marges indiquées ; zone grisée au-delà de fs/2, hors du domaine de validité du modèle moyen.
- `poles_boucle_fermee_zn.pdf` : pôles de la boucle fermée ZN aux 9 coins (disque unité, et agrandissement autour de z = 1 ; le pôle réel vers 0,72 est proche du pôle du filtre de dérivée, z = 1 − N·Tc = 0,709).
- `elm_marge_phase_ensemble.pdf` : pire marge de phase sur les 9 coins dans le plan (xP, xI), coupes xD = 0,5 ; 1 ; 2, frontière g = 0 en noir, départ (1 ; 1 ; 1) en croix sur la coupe xD = 1 ; en blanc, réglages non admissibles.
- `elm_distributions_ensemble.pdf` : distributions du pire module des pôles et de la pire marge de gain sur l'ensemble (grille 2^(1/16)), repère ZN.

