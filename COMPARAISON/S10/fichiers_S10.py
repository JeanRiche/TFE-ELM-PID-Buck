# =============================================================================
# fichiers_S10.py
#
# OBJECTIF
#   Ecrire l'essai S10 dans les quatre dossiers de methode (PSO_PID,
#   FUZZY_PID, ELM_PID, PINN_PID) sous la meme forme que les onze essais de
#   developpement, pour que le banc Python et Simulink le lisent comme S1 a
#   S9 :
#     - scenario_S10.mat, ecrit par les fonctions et l'appel savemat de
#       ../scenarios/scenarios_communs.py, recopies a l'identique
#       ci-dessous ;
#     - scenarios_communs.json : la liste "essais" (les onze essais de
#       developpement) ne change pas ; S10 est ajoute dans une liste a part,
#       "essais_complementaires", lue par banc_commun.py et par les scripts
#       Simuler_*.m. S10 n'entre pas dans le cout J.
#
# CONTROLES (le script s'arrete si l'un echoue)
#   1. Les fonctions recopiees redonnent, a l'identique, les profils de S3,
#      S4 et S5 deja ecrits dans ELM_PID (charge, Vin, consigne, bruit).
#   2. Les profils de S10 sont identiques, bit a bit, a ceux de
#      construire_S10() dans scenario_S10.py (celui qui a produit
#      resultats_S10.json) ; R0, q, duree et evenements aussi.
#   3. Relu par la fonction charger_scenario de banc_commun.py, le fichier
#      ecrit redonne les memes profils.
#
# PARTICULARITE DE S10
#   scenarios_communs.py refuse une conductance gx negative (la charge
#   electronique n'y fournit jamais de courant ; R0 est la plus grande
#   resistance de l'essai). S10 est defini avec R0 = 5 ohms et
#   gx = 1/R - 1/5 < 0 pour R = 20, 25 et 30 ohms. La source de courant
#   commandee de Buck_Commun.slx (i = gx x Vout, sans saturation) fournit
#   alors du courant, et la charge equivalente vaut bien 20, 25 ou 30 ohms.
#   Le controle gx >= 0 de scenarios_communs.py n'est donc pas applique a S10.
#
# COMMENT LANCER CE SCRIPT
#   python fichiers_S10.py      (quelques secondes)
# =============================================================================

import os
import json
import numpy as np
from scipy.io import savemat, loadmat

DOSSIER = os.path.dirname(os.path.abspath(__file__))
RACINE = os.path.dirname(os.path.dirname(DOSSIER))
DOSSIERS = ["PSO_PID", "FUZZY_PID", "ELM_PID", "PINN_PID"]


# --- Recopie a l'identique de scenarios/scenarios_communs.py, etapes 0 et 1 ---

TE = 1.0 / 220000.0                           # pas commun des profils = periode du regulateur Tc
VREF = 100.0                                  # consigne nominale (V)
VIN_NOM = 200.0                               # tension d'entree nominale (V)
R_NOM = 5.0                                   # charge nominale (ohms)
R_F2 = 5.0 + (-2.0 / 15.0) * 5.0              # charge de F2 (article) : 4.333 ohms
SANS_QUANTIF = 1e-9                           # pas "nul" : la quantification ne change rien


def pas_de(instant):
    """Numero du pas correspondant a un instant (en secondes)."""
    return int(round(instant / TE))


def palier(n, valeur_base, paliers):
    """Profil constant par morceaux : valeur_base partout, sauf sur chaque
    fenetre (debut, fin, valeur) de la liste 'paliers'. Les fenetres sont
    definies en numeros de pas entiers, pour que Simulink et le banc
    basculent exactement au meme pas."""
    x = np.full(n, float(valeur_base))
    for debut, fin, valeur in paliers:
        x[pas_de(debut):pas_de(fin)] = valeur
    return x


def conductance(r_equivalente, r0):
    """Conductance gx que doit tirer la charge electronique pour que la
    charge totale (R0 en parallele) vaille r_equivalente."""
    return 1.0 / r_equivalente - 1.0 / r0


def essai(code, nom, duree, R0, evenements, description, dvref=None, vin=None, gx=None, bruit=None, q=SANS_QUANTIF):
    """Rassemble un essai. Les profils absents prennent leur valeur nominale."""
    n = pas_de(duree) + 1                                  # nombre de pas, instant 0 compris
    return {"code": code, "nom": nom, "duree": duree, "R0": R0, "q": q,
            "evenements": list(evenements), "description": description,
            "t": np.arange(n) * TE,
            "dvref": np.zeros(n) if dvref is None else dvref(n),
            "vin": np.full(n, VIN_NOM) if vin is None else vin(n),
            "gx": np.zeros(n) if gx is None else gx(n),
            "bruit": np.zeros(n) if bruit is None else bruit(n)}


def vin_pm20(n):
    """Vin : 160 V sur [0.05 ; 0.07), 240 V sur [0.10 ; 0.12), 200 V ailleurs."""
    return palier(n, VIN_NOM, [(0.05, 0.07, 160.0), (0.10, 0.12, 240.0)])


def ecrire(e, dossier):
    """Meme appel que l'etape 3 de scenarios_communs.py."""
    savemat(os.path.join(dossier, f"scenario_{e['code']}.mat"),
            {"code": e["code"], "nom": e["nom"], "duree": e["duree"], "Te": TE, "R0": e["R0"], "q": e["q"],
             "evenements": np.array(e["evenements"], dtype=float).reshape(1, -1),
             "description": e["description"], "t": e["t"].reshape(-1, 1),
             "dvref": e["dvref"].reshape(-1, 1), "vin": e["vin"].reshape(-1, 1),
             "gx": e["gx"].reshape(-1, 1), "bruit": e["bruit"].reshape(-1, 1)},
            do_compression=True)


def lire(chemin):
    m = loadmat(chemin)
    return {"code": str(m["code"][0]), "nom": str(m["nom"][0]), "description": str(m["description"][0]),
            "duree": float(m["duree"].item()), "Te": float(m["Te"].item()), "R0": float(m["R0"].item()),
            "q": float(m["q"].item()), "evenements": [float(x) for x in np.ravel(m["evenements"])],
            "t": m["t"].ravel(), "dvref": m["dvref"].ravel(), "vin": m["vin"].ravel(), "gx": m["gx"].ravel(),
            "bruit": m["bruit"].ravel()}


PROFILS = ("dvref", "vin", "gx", "bruit")


def identiques(a, b, avec_t=True):
    ok = all(np.array_equal(a[k], b[k]) for k in PROFILS)
    ok = ok and a["R0"] == b["R0"] and a["q"] == b["q"] and a["duree"] == b["duree"]
    ok = ok and list(a["evenements"]) == list(b["evenements"])
    return ok and (np.array_equal(a["t"], b["t"]) if avec_t else len(a["t"]) == len(b["t"]))


# --- Controle 1 : les fonctions recopiees redonnent S3, S4 et S5 ---

TEMOINS = [
    essai("S3", "F2 (article)", 0.2, R_NOM, [0.05, 0.07], "charge de 4.333 ohms sur [0.05 ; 0.07) s",
          gx=lambda n: palier(n, 0.0, [(0.05, 0.07, conductance(R_F2, R_NOM))])),
    essai("S4", "Charge +-20 %", 0.2, 6.0, [0.05, 0.07, 0.10, 0.12],
          "4 ohms sur [0.05 ; 0.07) s, 6 ohms sur [0.10 ; 0.12) s, 5 ohms ailleurs",
          gx=lambda n: palier(n, conductance(R_NOM, 6.0), [(0.05, 0.07, conductance(4.0, 6.0)), (0.10, 0.12, 0.0)])),
    essai("S5", "Vin +-20 %", 0.2, R_NOM, [0.05, 0.07, 0.10, 0.12],
          "Vin = 160 V sur [0.05 ; 0.07) s, 240 V sur [0.10 ; 0.12) s", vin=vin_pm20)]
for e in TEMOINS:
    ref = lire(os.path.join(RACINE, "ELM_PID", f"scenario_{e['code']}.mat"))
    if not (identiques(e, ref) and ref["nom"] == e["nom"] and ref["description"] == e["description"]
            and ref["Te"] == TE):
        raise RuntimeError(f"Les fonctions recopiees ne redonnent pas {e['code']}.")
print("Controle 1 : les fonctions recopiees de scenarios_communs.py redonnent S3, S4 et S5 a l'identique.")


# --- S10 ---

R0_S10 = 5.0                                  # comme scenario_S10.py (gx < 0 apres 50 ms, voir l'en-tete)
S10 = essai("S10", "Changement durable du point", 0.2, R0_S10, [0.05, 0.10, 0.115, 0.13, 0.145, 0.16, 0.175],
            "essai complementaire, hors du cout J : Vin 200 -> 160 V et charge 5 -> 25 ohms a 50 ms (durable) ; "
            "20 ohms sur [0.10 ; 0.115) s, 30 ohms sur [0.13 ; 0.145) s, 20 ohms sur [0.16 ; 0.175) s",
            vin=lambda n: palier(n, VIN_NOM, [(0.05, 10.0, 160.0)]),
            gx=lambda n: palier(n, 0.0, [(0.05, 10.0, conductance(25.0, R0_S10)),
                                         (0.10, 0.115, conductance(20.0, R0_S10)),
                                         (0.13, 0.145, conductance(30.0, R0_S10)),
                                         (0.16, 0.175, conductance(20.0, R0_S10))]))

# --- Controle 2 : identique a construire_S10() de scenario_S10.py ---

with open(os.path.join(DOSSIER, "scenario_S10.py"), encoding="utf-8") as f:
    src = f.read()
src_fonction = src[src.index("def construire_S10():"):src.index("\nS10 = construire_S10()")]
H = 1.0 / (22000.0 * 1200.0)                  # constantes de banc_commun.py
ns = {"np": np, "TC": 120 * H}
exec(src_fonction, ns)
REF10 = ns["construire_S10"]()
if not identiques(S10, REF10, avec_t=False):
    raise RuntimeError("S10 ecrit ici differe de construire_S10() de scenario_S10.py.")
ecart_t = float(np.max(np.abs(S10["t"] - REF10["t"])))
print(f"Controle 2 : S10 identique, bit a bit, a construire_S10() de scenario_S10.py (Vin, gx, consigne, bruit, "
      f"R0, q, duree, evenements) ; instants t : ecart {ecart_t:.1e} s (t ne sert qu'aux figures et au "
      f"decalage d'un demi-pas de charger_scenario.m, calcule depuis Te).")
r_eq = 1.0 / (1.0 / S10["R0"] + S10["gx"])
for tms in (10, 60, 105, 120, 135, 150, 165, 190):
    k = pas_de(tms * 1e-3)
    print(f"  t = {tms:3d} ms : Vin = {S10['vin'][k]:5.1f} V, charge = {r_eq[k]:5.2f} ohms, gx = {S10['gx'][k]:+.4f} S")


# --- Ecriture dans les quatre dossiers et controle 3 ---

ENTREE = {k: S10[k] for k in ("code", "nom", "duree", "R0", "q", "evenements", "description")}
NOTE = ("Essais complementaires : simules et compares a Simulink comme les onze essais, mais hors du cout J "
        "(J reste defini sur les onze essais de la liste 'essais').")
for d in DOSSIERS:
    dossier = os.path.join(RACINE, d)
    ecrire(S10, dossier)
    relu = lire(os.path.join(dossier, "scenario_S10.mat"))
    if not (identiques(S10, relu) and relu["Te"] == TE and relu["code"] == "S10"):
        raise RuntimeError(f"{d} : scenario_S10.mat relu ne redonne pas S10.")
    chemin = os.path.join(dossier, "scenarios_communs.json")
    with open(chemin, encoding="utf-8") as f:
        texte = f.read()
    liste = json.loads(texte)
    if [e["code"] for e in liste["essais"]] != ["S1", "S2", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]:
        raise RuntimeError(f"{d} : scenarios_communs.json ne contient pas les onze essais attendus.")
    base = {k: v for k, v in liste.items() if k not in ("essais_complementaires", "note_essais_complementaires")}
    liste = dict(base, note_essais_complementaires=NOTE, essais_complementaires=[ENTREE])
    nouveau = json.dumps(liste, indent=1)
    if not nouveau.startswith(json.dumps(base, indent=1)[:-2]):
        raise RuntimeError(f"{d} : la partie existante de scenarios_communs.json changerait.")
    with open(chemin, "w", encoding="utf-8") as f:
        f.write(nouveau)
    print(f"{d} : scenario_S10.mat ecrit et relu a l'identique ; scenarios_communs.json : S10 ajoute aux essais "
          f"complementaires (les onze essais inchanges).")
