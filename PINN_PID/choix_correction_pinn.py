# =============================================================================
# choix_correction_pinn.py
#
# OBJECTIF
#   Choisir mecaniquement, sur les seuls essais de reglage E1 a E4 (points R,
#   Vin absents des onze essais de developpement), la combinaison de
#   corrections du PINN-PID, selon la regle de regle_choix_correction_pinn.json :
#     1. controles : gradient de l'horizon (c = 0 et c = 1, faits au
#        chargement du banc), fonctions et table de l'ensemble admissible
#        identiques a celles de l'ELM-PID, identite sans adaptation (ecart 0
#        sur d au PID de Ziegler-Nichols sur E1 a E4) ;
#     2. les 8 combinaisons (aucune, C1, C2, C3, C1+C2, C1+C3, C2+C3,
#        C1+C2+C3) sur E1 a E4 : IAE du demarrage et apres 30 ms, rapports a
#        Ziegler-Nichols, critere, retours dans +-1 V, gains ;
#     3. application de la regle ; ecriture de choix_correction_pinn.json,
#        lu ensuite par banc_pinn_pid.py.
#   Rien n'est execute sur les onze essais de developpement ici.
#
# COMMENT LANCER CE SCRIPT
#   python choix_correction_pinn.py      (quelques minutes)
# =============================================================================

import os
import json
import hashlib
import numpy as np

try:
    DOSSIER_PINN = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PINN = os.getcwd()
__file__ = os.path.join(DOSSIER_PINN, "banc_pinn_pid.py")
with open(__file__, encoding="utf-8") as f:
    _SRC = f.read()
exec(_SRC[:_SRC.index("# %% ETAPE 2 :")])   # banc commun, PINN, observateur, horizon, PINNPID, controles

with open(os.path.join(DOSSIER_PINN, "regle_choix_correction_pinn.json"), "rb") as f:
    _octets = f.read()
REGLE = json.loads(_octets.decode("utf-8"))
print(f"Regle de choix lue : regle_choix_correction_pinn.json, sha256 {hashlib.sha256(_octets).hexdigest()}")


# %% ETAPE 1 : essais de reglage E1 a E4 (memes definitions que dans recherche_pso_pid.py)

def essai_reglage(code, R0, vin0, duree, evenements=(), gx=None, vin=None, dvref=None, bruit=False, q=1e-9):
    n = int(round(duree / TC)) + 1
    t = np.arange(n) * TC

    def palier(base, morceaux):
        y = np.full(n, base, dtype=float)
        for t0, t1, val in morceaux:
            y[(t >= t0 - 1e-12) & (t < t1 - 1e-12)] = val
        return y

    sc = {"code": code, "nom": code, "duree": duree, "R0": R0, "q": q, "evenements": list(evenements), "t": t,
          "dvref": palier(0.0, dvref or []), "vin": palier(vin0, vin or []), "gx": palier(0.0, gx or []),
          "bruit": np.zeros(n)}
    if bruit:                                  # meme nature que S7b (10 mV, 10 kHz), autre graine (7)
        blanc = np.random.default_rng(7).standard_normal(n)
        a = np.exp(-2.0 * np.pi * 10e3 * TC)
        y, etat = np.zeros(n), 0.0
        for k in range(n):
            etat = a * etat + (1.0 - a) * blanc[k]
            y[k] = etat
        sc["bruit"] = 0.01 * y / np.std(y)
    return sc


E = [essai_reglage("E1", 7.0, 190.0, 0.090, (0.040, 0.065), gx=[(0.040, 0.065, 1 / 4.5 - 1 / 7.0)]),
     essai_reglage("E2", 15.0, 210.0, 0.090, (0.040, 0.065), vin=[(0.040, 0.065, 175.0), (0.065, 1.0, 235.0)]),
     essai_reglage("E3", 60.0, 220.0, 0.090, (0.040, 0.065), dvref=[(0.040, 0.065, 10.0)]),
     essai_reglage("E4", 9.0, 200.0, 0.060, (), bruit=True, q=2.5 / 4096 * 80)]
for sc in E:                                   # aucun point de reglage n'est un point des onze essais
    for s_ in SCENARIOS:
        if abs(s_["R0"] - sc["R0"]) < 1e-9 and np.any(np.abs(s_["vin"] - sc["vin"][0]) < 1e-9):
            raise RuntimeError(f"Le point {sc['code']} apparait dans l'essai {s_['code']}.")
ESSAIS_E = {sc["code"]: sc for sc in E}

titre("ETAPE 1 : identite sans adaptation (E1 a E4)")
zn = evaluer(lambda: PIDClassique(), essais=ESSAIS_E)
ECART_ID = 0.0
for corr in ((), ("C1", "C2", "C3")):
    fige = evaluer(lambda: PINNPID(adapter=False, corrections=corr), essais=ESSAIS_E)
    ECART_ID = max(ECART_ID, max(float(np.max(np.abs(fige[c]["sim"]["d"] - zn[c]["sim"]["d"]))) for c in ESSAIS_E))
print(f"  PINN-PID sans adaptation (sans correction et avec C1+C2+C3) contre Ziegler-Nichols, ecart maximal "
      f"sur d : {ECART_ID:.1e}")
if ECART_ID != 0.0:
    raise RuntimeError("Sans adaptation, le PINN-PID ne redonne pas le PID de Ziegler-Nichols.")


# %% ETAPE 2 : les 8 combinaisons sur E1 a E4

titre("ETAPE 2 : les 8 combinaisons sur E1 a E4 (IAE en mV.s)")


def iae(r):
    e = np.abs(r["sim"]["consigne"] - r["sim"]["v"])
    return float(np.sum(e[:K30]) * TC), float(np.sum(e[K30:]) * TC)


def nom_comb(c):
    return "+".join(c) if c else "aucune"


T_ZN = {c: iae(zn[c]) for c in ESSAIS_E}
TABLEAU = {}
for comb in REGLE["combinaisons"]:
    comb = tuple(comb)
    res = evaluer(lambda: PINNPID(corrections=comb), essais=ESSAIS_E)
    rapports, lignes, non, hors = [], [], [], 0
    for sc in E:
        c = sc["code"]
        dem, apres = iae(res[c])
        rapports += [apres / T_ZN[c][1], dem / T_ZN[c][0]]
        if not res[c]["revenus"]:
            non.append(c)
        jr = res[c]["reg"].journal
        X = np.array([x["x"] for x in jr])
        hors += sum(1 for x in jr if x["admis"] is False)
        n_opt = sum(1 for x in jr if x["adapte"])
        lignes.append(f"    {c} demarrage {dem * 1e3:6.2f} ({dem / T_ZN[c][0]:.3f}), apres 30 ms {apres * 1e3:6.2f} "
                      f"({apres / T_ZN[c][1]:.3f}) ; optimise sur {n_opt:3d}/{len(jr)} fenetres ; x final "
                      f"({X[-1, 0]:.3f}, {X[-1, 1]:.3f}, {X[-1, 2]:.3f}) ; apres 30 ms P {X[59:, 0].min():.2f}-"
                      f"{X[59:, 0].max():.2f}, I {X[59:, 1].min():.2f}-{X[59:, 1].max():.2f}, D "
                      f"{X[59:, 2].min():.2f}-{X[59:, 2].max():.2f}")
    crit = float(np.mean(rapports))
    TABLEAU[nom_comb(comb)] = {"corrections": list(comb), "critere": crit,
                               "apres_30ms": float(np.mean(rapports[0::2])), "demarrage": float(np.mean(rapports[1::2])),
                               "rapports": rapports, "non_revenus": non, "fenetres_hors_ensemble": hors}
    print(f"  {nom_comb(comb):9s} critere {crit:.4f} (apres 30 ms {np.mean(rapports[0::2]):.3f}, demarrage "
          f"{np.mean(rapports[1::2]):.3f}) ; non revenus : {','.join(non) or '-'}"
          + (f" ; fenetres hors de l'ensemble : {hors}" if "C3" in comb else ""), flush=True)
    for l_ in lignes:
        print(l_, flush=True)


# %% ETAPE 3 : application mecanique de la regle

titre("ETAPE 3 : choix (regle_choix_correction_pinn.json)")
admis = {n: t for n, t in TABLEAU.items() if not t["non_revenus"]}
meilleur = min(t["critere"] for t in admis.values())
proches = {n: t for n, t in admis.items() if t["critere"] <= (1.0 + REGLE["tolerance_relative"]) * meilleur}
retenue = min(proches, key=lambda n: (len(proches[n]["corrections"]), proches[n]["critere"]))
print(f"  Meilleur critere {meilleur:.4f} ; a 1 % pres : " + ", ".join(f"{n} ({t['critere']:.4f})" for n, t in proches.items()))
print(f"  Combinaison retenue : {retenue} (critere {TABLEAU[retenue]['critere']:.4f})")
with open(os.path.join(DOSSIER_PINN, "choix_correction_pinn.json"), "w", encoding="utf-8") as f:
    json.dump({"regle_sha256": hashlib.sha256(_octets).hexdigest(),
               "combinaison_retenue": TABLEAU[retenue]["corrections"], "nom": retenue,
               "meilleur_critere": meilleur, "a_1_pct": list(proches), "ecart_identite_ZN": ECART_ID,
               "ecart_gradient": ECART_GRAD, "ecart_gradient_saturation": ECART_GRAD_SAT,
               "ecart_gradient_C2": ECART_GRAD_C2, "tableau": TABLEAU}, f, indent=1)
print("  choix_correction_pinn.json ecrit.")
