# Annexe du mémoire : régulation d'un convertisseur Buck par un PID auto-adaptatif à ELM

Ce dépôt contient les scripts, les données et les résultats cités dans le mémoire
« Optimisation de la régulation d'un convertisseur Buck par un contrôleur PID
auto-adaptatif basé sur l'Extreme Learning Machine (ELM) ». Il permet de recalculer les
chiffres, les tableaux et les figures des chapitres 3, 4 et 5, qui comparent cinq
régulateurs du même convertisseur : PID de Ziegler-Nichols (ZN), PSO-PID, Fuzzy-PID,
ELM-PID et PINN-PID.

Deux outils de calcul se complètent.

- Le banc Python simule le convertisseur et les régulateurs au pas du régulateur
  (T_c = 1/220 000 s). Il produit les résultats chiffrés (IAE, métriques, séparation,
  ablation, stabilité).
- Les modèles Simulink (MATLAB R2024a) simulent le même circuit avec les blocs
  Simscape Electrical. Ils produisent les figures du chapitre 5 et la validation croisée
  du chapitre 3. Les scripts MATLAB et les modèles `.slx` ne peuvent pas être exécutés
  sans MATLAB.

## Arborescence

| Dossier | Rôle |
|---|---|
| `ELM_PID/` | ELM-PID (bloc PID commun à gains externes) : banc, apprentissage hors ligne, ensemble admissible, adaptateur MATLAB, construction et simulation Simulink, données d'apprentissage. Contient aussi le banc commun `banc_commun.py` et les essais (`scenario_*.mat`, `scenarios_communs.json`) |
| `ELM_PID_INCREMENTAL/` | ELM-PID à loi incrémentale de Lu : banc Python et figure de comparaison avec le bloc PID commun (chapitre 4) |
| `PSO_PID/` | PSO-PID : recherche par essaim (`recherche_pso_pid.py`), banc, construction et simulation Simulink. Contient aussi `Simuler_ZN.m` (Ziegler-Nichols sous Simulink) |
| `FUZZY_PID/` | Fuzzy-PID (ordonnancement flou de Zhao, Tomizuka et Isaka) : banc, bloc MATLAB, construction et simulation Simulink |
| `PINN_PID/` | PINN-PID (optimisation d'Ito et Wasa en temps réel) : apprentissage du réseau, boîte de gains, choix de la correction, banc, bloc MATLAB, construction et simulation Simulink |
| `COMPARAISON/` | Comparaison des cinq méthodes : métriques, scénario S10, séparation des méthodes, ablation du réseau, classement par essai, comparaison des deux lois de l'ELM-PID, figures du chapitre 5, génération des essais |
| `STABILITE/` | Analyse de stabilité linéaire à gains figés (tableaux et figures des chapitres 3 et 4) |
| `REDACTION/` | Schémas des chapitres 3 et 4, relevé de Ziegler-Nichols, recalcul de quatre valeurs du mémoire |

Chaque dossier de méthode contient ses propres copies des fichiers d'essais
(`scenario_*.mat`, `scenarios_communs.json`) et du banc commun, car les scripts les lisent
dans leur dossier. Ces copies sont identiques entre dossiers ; seul `ELM_PID_INCREMENTAL/`
ne contient pas S10.

Types de fichiers :

- `banc_*.py`, `recherche_*.py`, `entrainement_*.py` : calculs Python ;
- `*_resultats.json`, `predictions_banc_*.json` : résultats du banc, lus par les autres
  scripts et par les scripts MATLAB (qui comparent Simulink au banc) ;
- `*.mat` : essais, modèles appris, réglages et données d'apprentissage ;
- `*.m` et `Buck_Commun.slx` : MATLAB et Simulink (modèle commun du convertisseur).

## Prérequis

Python 3.10 ou plus (résultats vérifiés avec Python 3.13, NumPy 2.5, SciPy 1.18,
Matplotlib 3.11, schemdraw 0.23) :

```
pip install -r requirements.txt
```

Un système LaTeX (`latex`, `dvipng`, paquets `newtx` et `icomma`) est nécessaire pour les
figures qui composent leur texte avec LaTeX : `ELM_PID_INCREMENTAL/figure_loi_incrementale_bloc_commun.py`,
`REDACTION/chapitre3/figures/schema_chapitre3.py` et `REDACTION/chapitre4/figures/schemas_chapitre4.py`.
Le paquet `schemdraw` ne sert qu'à `schema_chapitre3.py`.

MATLAB R2024a avec Simulink, Simscape et Simscape Electrical (Specialized Power Systems)
pour les scripts `.m` et le modèle `Buck_Commun.slx`. Les scripts se lancent depuis le
dossier de la méthode (dossier courant de MATLAB), sauf `Metriques_Simulink.m` et
`Figures_Chapitre5.m`, qui se placent dans leur propre dossier.

Durées indicatives (processeur à quatre cœurs) : banc ELM-PID 1,5 min, PSO-PID 1 min,
Fuzzy-PID 4 min, PINN-PID 10 min ; `metriques_banc.py` quelques minutes ;
`stabilite.py` 30 à 40 min ; `recherche_pso_pid.py` avec `RECALCULER = True` 35 à 70 min
(par défaut, il reprend `recherche_pso_pid.json`).

## Reproduction : ordre d'exécution

Les commandes se lancent depuis le dossier indiqué. Les fichiers déjà présents
(`*.mat`, `*.json`) sont les résultats de ces calculs ; relancer un script les remplace.

1. Essais communs : `python scenarios_communs.py` dans `COMPARAISON/scenarios/` produit
   `scenario_S1.mat` à `scenario_S9.mat` et `scenarios_communs.json`, identiques à ceux
   des dossiers de méthode. `python fichiers_S10.py` dans `COMPARAISON/S10/` ajoute S10 à
   ces dossiers.
2. ELM-PID, dans `ELM_PID/` : `excitation_donnees_elm.py`, puis (données d'apprentissage
   Simulink : `Construction_Modele_Donnees_ELM.m` et `Generer_Donnees_ELM.m`)
   `entrainement_elm.py`, `ensemble_gains_elm.py`, `banc_elm_pid.py`.
3. PSO-PID, dans `PSO_PID/` : `recherche_pso_pid.py`, `banc_pso_pid.py`.
4. Fuzzy-PID, dans `FUZZY_PID/` : `banc_fuzzy_pid.py`.
5. PINN-PID, dans `PINN_PID/` : `entrainement_pinn.py`, `boite_gains_pinn.py`,
   `choix_correction_pinn.py`, puis `banc_pinn_pid.py`. Copier avant
   `ELM_PID/ensemble_gains_elm.mat` dans `PINN_PID/` (le banc vérifie que les deux
   fichiers sont identiques).
6. Copier `banc_elm_pid_resultats.json`, `banc_pso_pid_resultats.json`,
   `banc_fuzzy_pid_resultats.json` et `banc_pinn_pid_resultats.json` dans `COMPARAISON/`.
7. Comparaison, dans `COMPARAISON/` : `S10/scenario_S10.py`, `metriques/metriques_banc.py`,
   `separation/separation_methodes.py`, `ablation_S10/ablation_reseau_elm.py`,
   `comparaison_methodes.py`.
8. Simulink, dans chaque dossier de méthode : `Construction_<méthode>.m`, puis
   `Simuler_<méthode>.m` (validation sur les onze essais), puis
   `Construire_Modeles_Comparaison.m` et `Simuler_Modeles_Comparaison.m` (S1, S2, S3, S8a,
   S10). Pour Ziegler-Nichols : `PSO_PID/Simuler_ZN.m`, et
   `Construire_Modeles_Comparaison.m` avec `AVEC_ZIEGLER_NICHOLS = true`.

## Correspondance avec le mémoire

Les numéros de figures et de tableaux sont ceux du mémoire.

| Élément du mémoire | Commande | Fichier produit |
|---|---|---|
| Fig. 3.1, schéma du convertisseur | `python REDACTION/chapitre3/figures/schema_chapitre3.py` | `REDACTION/chapitre3/figures/Fig3_1_convertisseur_Buck.pdf` |
| Fig. 3.3, modèle Simulink commun | `Exporter_Buck_Commun_Times.m` (MATLAB, dossier `ELM_PID/`) | `ELM_PID/Buck_Commun.pdf` |
| Tab. 3.2 à 3.4, dimensionnement : gain statique avec pertes, ESR tolérée, limites de la conduction continue ; pic de l'action dérivée (chap. 4) | `python REDACTION/verifications_chiffres_memoire.py` | affichage console (code de sortie 0 si les valeurs coïncident avec le texte) |
| Fig. 3.4, relevé de Ziegler-Nichols ; gains P, I, D, N | `Figure_Releve_Ziegler_Nichols.m` (MATLAB, dossier `REDACTION/chapitre3/figures/`) | `figures_zn/fig3_02_releve_Ziegler_Nichols.*` |
| Tab. 3.5 à 3.7, fig. 3.5 et 3.6 : pôles, boucle fermée de ZN aux neuf coins, charge légère | `python STABILITE/stabilite.py` | `STABILITE/resultats_stabilite.md`, `.json`, `figures/bode_boucle_ouverte_zn_nominal.pdf`, `figures/poles_boucle_fermee_zn.pdf` |
| Tab. 3.8 et 3.9, essais de développement S1 à S9 et essais de réglage E1 à E4 | `python COMPARAISON/scenarios/scenarios_communs.py` ; `PSO_PID/recherche_pso_pid.py` (définition de E1 à E4) | `scenario_*.mat`, `scenarios_communs.json` ; `PSO_PID/recherche_pso_pid.json` |
| Tab. 3.10, écart d'IAE entre Simulink et le banc | `Simuler_ZN.m`, `Simuler_PSO_PID.m`, `Simuler_Fuzzy_PID.m`, `Simuler_ELM_PID.m`, `Simuler_PINN_PID.m` (MATLAB) | `resultats_Buck_Commun_<méthode>.mat` ; comparaison affichée aux `predictions_banc_*.json` |
| Fig. 4.1 à 4.3, schémas de l'ELM-PID | `python REDACTION/chapitre4/figures/schemas_chapitre4.py` | `REDACTION/chapitre4/figures/Fig4_1_*.pdf` à `Fig4_3_*.pdf` |
| Apprentissage hors ligne (chap. 4) : ordre du modèle (rapports 1,193 et 1,078), choix du réseau, jacobien vrai par classe de charge | `python ELM_PID/entrainement_elm.py` | `ELM_PID/elm_pid_modele.mat`, `ELM_PID/entrainement_elm_resultats.json` |
| Ensemble admissible (marge de départ, coupure, table) | `python ELM_PID/ensemble_gains_elm.py` | `ELM_PID/ensemble_gains_elm.mat`, `.json` |
| Tab. 4.4 et fig. 4.4, loi incrémentale de Lu contre bloc PID commun | `python COMPARAISON/loi_incrementale_contre_bloc_commun.py` ; `python ELM_PID_INCREMENTAL/figure_loi_incrementale_bloc_commun.py` | `COMPARAISON/loi_incrementale_contre_bloc_commun.json` ; `ELM_PID_INCREMENTAL/figures/fig4_05_loi_incrementale_bloc_commun.pdf` |
| Fig. 4.5, pire marge de phase sur l'ensemble admissible | `python STABILITE/stabilite.py` | `STABILITE/figures/elm_marge_phase_ensemble.pdf` |
| Fig. 4.6, modèle Simulink de l'ELM-PID | `Construction_ELM_PID.m`, puis `Exporter_ELM_PID_Times.m` (MATLAB, dossier `ELM_PID/`) | `ELM_PID/fig4_06_modele_Simulink_ELM_PID.pdf` |
| PINN-PID (chap. 4 et 5) : réseau, boîte de gains, correction retenue | `python PINN_PID/entrainement_pinn.py`, `boite_gains_pinn.py`, `choix_correction_pinn.py` | `PINN_PID/pinn_pid_modele.mat`, `boite_gains_pinn.json`, `choix_correction_pinn.json` |
| IAE, dépassements, établissements, évènements de chaque méthode (chap. 5) | `python ELM_PID/banc_elm_pid.py`, `PSO_PID/banc_pso_pid.py`, `FUZZY_PID/banc_fuzzy_pid.py`, `PINN_PID/banc_pinn_pid.py` | `banc_<méthode>_pid_resultats.json` dans chaque dossier |
| Tab. 5.1 et 5.2, paramètres et scénarios | `banc_commun.py` (constantes) ; `scenarios_communs.py`, `fichiers_S10.py` | `scenarios_communs.json` |
| Tab. 5.4 (ZN et ELM-PID), tab. 5.7 (temps de calcul par appel), mesures citées dans le texte | `python COMPARAISON/metriques/metriques_banc.py` | `COMPARAISON/metriques/metriques_banc.md`, `metriques_banc.csv` |
| Tab. 5.3 et 5.5, scénario S10 : gains figés à 50 ms, gains à 100 ms | `python COMPARAISON/S10/scenario_S10.py` | `COMPARAISON/S10/resultats_S10.json` |
| Tab. 5.6, séparation des méthodes (groupes ordonnés, IAE) | `python COMPARAISON/separation/separation_methodes.py` | `COMPARAISON/separation/separation_resultats.json` |
| Ablation du réseau (jacobien constant, sans OS-ELM ; tab. 5.8 et 5.9) | `python COMPARAISON/ablation_S10/ablation_reseau_elm.py` | `COMPARAISON/ablation_S10/ablation_resultats.json` |
| Classement par essai (annexe) | `python COMPARAISON/comparaison_methodes.py` | `COMPARAISON/comparaison_resultats.json` |
| Fig. 5.1 à 5.15 | `Figures_Chapitre5.m` (MATLAB, dossier `COMPARAISON/figures_chapitre5/`), après l'étape 8 | `COMPARAISON/figures_chapitre5/pdf_memoire/*.pdf` et `Fig5_*.png` |
| Écarts d'IAE Simulink - banc sur S1, S2, S3, S8a, S10 (chap. 5) | `Figures_Chapitre5.m` (contrôle affiché), `COMPARAISON/metriques/Metriques_Simulink.m` | console ; `metriques_simulink_comparaison.csv` et `.md` |

## Contrôles intégrés

- Chaque banc recalcule le PID de Ziegler-Nichols à l'identique et s'arrête si le régulateur
  sans adaptation s'en écarte.
- `metriques_banc.py`, `separation_methodes.py` et `ablation_reseau_elm.py` comparent leurs
  IAE nominales aux fichiers de résultats (écart relatif nul ou inférieur à 1e-9) et
  s'arrêtent sinon.
- `figure_loi_incrementale_bloc_commun.py` compare les valeurs du tableau 4.4 aux valeurs
  imprimées dans le mémoire avant de tracer.

## Remarques

- Les fichiers `.mat` et `.json` présents sont les résultats des scripts, recalculés à
  partir de ces mêmes scripts. Les temps de calcul (`duree_s`) varient d'une machine à
  l'autre.
- Les temps d'appel du régulateur de `metriques_banc.md` et `metriques_banc.csv`
  (tableau 5.7) sont ceux de la machine de calcul du mémoire. Sur une autre machine,
  toutes les autres valeurs sont retrouvées à l'identique, mais les temps changent ;
  seul leur rapport entre méthodes a un sens.
- Les modèles Simulink construits par les scripts (`Buck_Commun_<méthode>.slx`,
  `<méthode>_<scénario>.slx`) et leurs résultats (`resultats_*.mat`) ne sont pas suivis par git.
- L'ELM-PID à loi incrémentale n'a pas de modèle Simulink dans ce dépôt : sa comparaison
  avec le bloc PID commun repose sur le banc Python.
