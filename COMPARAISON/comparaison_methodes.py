# =============================================================================
# comparaison_methodes.py
#
# OBJECTIF
#   A partir des resultats du banc commun de chaque methode, classer les
#   cinq methodes sur chaque essai par l'IAE de la fenetre de classement
#   (30 ms a la fin) et ecrire le tableau par essai (comparaison_resultats.json)
#   et la figure d'annexe classement_par_essai.png.
#
# FICHIERS NECESSAIRES (dans ce dossier)
#   banc_pso_pid_resultats.json, banc_fuzzy_pid_resultats.json,
#   banc_elm_pid_resultats.json, banc_pinn_pid_resultats.json
#
# COMMENT LANCER CE SCRIPT
#   python comparaison_methodes.py      (quelques secondes)
# =============================================================================

import os
import json
import numpy as np

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()


def lire(nom):
    with open(os.path.join(DOSSIER, nom), encoding="utf-8") as f:
        return json.load(f)["references"]


SOURCES = {"PSO-PID": ("banc_pso_pid_resultats.json", "PSO-PID"),
           "Fuzzy-PID": ("banc_fuzzy_pid_resultats.json", "Fuzzy-PID"),
           "ELM-PID": ("banc_elm_pid_resultats.json", "ELM-PID"),
           "PINN-PID": ("banc_pinn_pid_resultats.json", "PINN-PID")}
REF = {n: lire(f) for n, (f, _) in SOURCES.items()}
DONNEES = {"Ziegler-Nichols": REF["PSO-PID"]["Ziegler-Nichols"]}
for n, (f, cle) in SOURCES.items():
    DONNEES[n] = REF[n][cle]
METHODES = list(DONNEES)                      # ZN, PSO, Fuzzy, ELM, PINN
CODES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]
CODES_IAE = ["S1", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]
CODES_DEM = ["S1", "S8a", "S8b"]


def titre(t):
    print("\n" + "=" * 100 + "\n " + t + "\n" + "=" * 100)


# %% ETAPE 1 : controle des donnees

titre("ETAPE 1 : controle des donnees")
ecart_zn = 0.0
for n in SOURCES:
    for c in CODES:
        a, b = REF[n]["Ziegler-Nichols"][c], DONNEES["Ziegler-Nichols"][c]
        ecart_zn = max(ecart_zn, abs(a["grandeurs"]["IAE"] - b["grandeurs"]["IAE"]) / b["grandeurs"]["IAE"],
                       abs(a["IAE_dem"] - b["IAE_dem"]) / b["IAE_dem"])
print(f"  Ziegler-Nichols identique dans les quatre fichiers : ecart relatif maximal {ecart_zn:.1e}")
if ecart_zn > 1e-9:
    raise RuntimeError("Les quatre bancs n'ont pas les memes grandeurs pour Ziegler-Nichols.")



# %% ETAPE 2 : chaque essai

def iae(m, c):
    return DONNEES[m][c]["grandeurs"]["IAE"]


def revenus(m, c):
    return DONNEES[m][c]["revenus"]


def ordre_essai(c):
    """Evenements revenus d'abord, puis IAE croissante."""
    return sorted(METHODES, key=lambda m: (not revenus(m, c), iae(m, c)))


titre("ETAPE 2 : IAE de 30 ms a la fin (mV.s) par essai et par methode ; classement")
print(f"  {'essai':5s}" + "".join(f"{m:>16s}" for m in METHODES) + "  classement")
TABLE = {}
for c in CODES:
    o = ordre_essai(c)
    TABLE[c] = {"ordre": o, "IAE_mVs": {m: iae(m, c) * 1e3 for m in METHODES},
                "revenus": {m: revenus(m, c) for m in METHODES}}
    ligne = f"  {c:5s}" + "".join(f"{iae(m, c) * 1e3:15.2f}{'' if revenus(m, c) else '*'}{' ' if revenus(m, c) else ''}"
                                   for m in METHODES)
    print(ligne + "  " + " > ".join(o))
print("  (* : un evenement au moins ne revient pas dans +-1 V ; classement : meilleur a gauche)")
print("\n  Demarrage (IAE de 0 a 30 ms, mV.s) :")
for c in CODES_DEM:
    print(f"    {c:4s}" + "".join(f"  {m} {DONNEES[m][c]['IAE_dem'] * 1e3:6.2f}" for m in METHODES))


with open(os.path.join(DOSSIER, "comparaison_resultats.json"), "w", encoding="utf-8") as f:
    json.dump({"essais": TABLE}, f, indent=1)
print("\n  comparaison_resultats.json ecrit.")


# %% ETAPE 5 : figure d'annexe (IAE rapportee a Ziegler-Nichols, essai par essai)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = np.array([[iae(m, c) / iae("Ziegler-Nichols", c) for c in CODES] for m in METHODES])
fig, ax = plt.subplots(figsize=(13, 3.8))
im = ax.imshow(np.log2(M), cmap="RdYlGn_r", vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(len(CODES)))
ax.set_xticklabels(CODES, fontsize=9)
ax.set_yticks(range(len(METHODES)))
ax.set_yticklabels(METHODES, fontsize=9)
for i, m in enumerate(METHODES):
    for j, c in enumerate(CODES):
        rang = TABLE[c]["ordre"].index(m) + 1
        ax.text(j, i, f"{M[i, j]:.2f}\n({rang})", ha="center", va="center", fontsize=7.5)
ax.set_title("IAE de 30 ms a la fin rapportee a Ziegler-Nichols (rang dans l'essai entre parentheses ; "
             "vert : meilleur que ZN)", fontsize=10, loc="left")
fig.colorbar(im, ax=ax, label="log2 du rapport")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER, "classement_par_essai.png"), dpi=110)
plt.close(fig)
print("  classement_par_essai.png ecrit.")
