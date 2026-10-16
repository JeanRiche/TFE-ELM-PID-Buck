# =============================================================================
# loi_incrementale_contre_bloc_commun.py
#
# OBJECTIF
#   Comparer, pour le chapitre 4, les deux formes de l'ELM-PID, chacune
#   prise dans son dossier :
#     - ../ELM_PID_INCREMENTAL : loi PID incrementale de Lu ;
#     - ../ELM_PID             : bloc PID commun a gains externes.
#   Seule la loi differe : modele ELM, OS-ELM, gradient, porte, zones
#   mortes, ensemble admissible (recalcule pour chaque loi) et projection
#   sont les memes, avec les memes reglages.
#   Essais : S4 pour la comparaison, avec la decomposition de J ; S2 et S1
#   pour montrer les deux mecanismes (a-coup de consigne, demarrage). On
#   ajoute, en reference, le PID de Ziegler-Nichols (point de depart des
#   deux formes).
#
# CE QUE PRODUIT CE SCRIPT
#   loi_incrementale_contre_bloc_commun.json et
#   loi_incrementale_contre_bloc_commun.png.
#
# COMMENT LANCER CE SCRIPT
#   python loi_incrementale_contre_bloc_commun.py
#   (une minute ; les dossiers ELM_PID et ELM_PID_INCREMENTAL a cote).
# =============================================================================

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
RACINE = os.path.dirname(DOSSIER)


def charger(dossier):
    """Definitions du banc d'une forme (jusqu'a l'etape 2), dans un espace
    de noms propre."""
    chemin = os.path.join(RACINE, dossier, "banc_elm_pid.py")
    ns = {"__file__": chemin, "__name__": "definitions"}
    with open(chemin, encoding="utf-8") as f:
        src = f.read()
    exec(src[:src.index("# %% ETAPE 2 :")], ns)
    with open(os.path.join(RACINE, dossier, "banc_elm_pid_resultats.json"), encoding="utf-8") as f:
        ns["RESULTATS"] = json.load(f)
    return ns


FORMES = {"incrementale (Lu)": charger("ELM_PID_INCREMENTAL"), "bloc PID commun": charger("ELM_PID")}
NS_B = FORMES["bloc PID commun"]
K30 = int(round(0.03 / NS_B["TC"]))


def lancer(ns, fabrique, code):
    sc = {s["code"]: s for s in ns["SCENARIOS"]}[code]
    reg = fabrique()
    sim = ns["simuler"](reg, sc)
    o = ns["grandeurs"](sim, sc)
    e = sim["consigne"] - sim["v"]
    return {"o": o, "sim": sim, "reg": reg, "IAE_dem": float(np.sum(np.abs(e[:K30])) * ns["TC"])}


# 1. Cout J et sa decomposition (fichiers de resultats de chaque dossier)
print("Cout J (13 termes) et decomposition, d'apres banc_elm_pid_resultats.json de chaque dossier")
TABLE_J = {}
for nom, ns in FORMES.items():
    r = ns["RESULTATS"]
    TABLE_J[nom] = {"J": r["J"]["ELM-PID"], "apres30": r["decomposition"]["ELM-PID"][0],
                    "demarrage": r["decomposition"]["ELM-PID"][1]}
    print(f"  {nom:22s} J {TABLE_J[nom]['J']:.3f} ; apres 30 ms {TABLE_J[nom]['apres30']:.3f} ; "
          f"demarrage {TABLE_J[nom]['demarrage']:.3f}")

# 2. Les trois essais
CODES = ("S4", "S2", "S1")
RES = {}
for code in CODES:
    RES[code] = {"Ziegler-Nichols": lancer(NS_B, lambda: NS_B["PIDClassique"](), code)}
    for nom, ns in FORMES.items():
        RES[code][nom] = lancer(ns, lambda ns=ns: ns["ELMPID"](), code)
    for nom, ns in FORMES.items():                       # controle : le banc redonne son fichier
        attendu = ns["RESULTATS"]["references"]["ELM-PID"][code]["grandeurs"]["IAE"]
        if abs(RES[code][nom]["o"]["IAE"] - attendu) > 1e-12:
            raise RuntimeError(f"{nom}, {code} : IAE {RES[code][nom]['o']['IAE']} au lieu de {attendu}.")

print("\nIAE (mV.s) et grandeurs ; controle : chaque forme redonne l'IAE de son fichier de resultats")
SORTIE = {"J": TABLE_J, "essais": {}}
for code in CODES:
    print(f"  {code}")
    SORTIE["essais"][code] = {}
    for nom, r in RES[code].items():
        o = r["o"]
        ev = o["evenements"]
        x_fin = (r["reg"].journal[-1]["x"].tolist() if hasattr(r["reg"], "journal") else [1.0, 1.0, 1.0])
        nb = (sum(1 for j in r["reg"].journal if j["change"]) if hasattr(r["reg"], "journal") else 0)
        SORTIE["essais"][code][nom] = {"IAE_apres30_mVs": o["IAE"] * 1e3, "IAE_dem_mVs": r["IAE_dem"] * 1e3,
                                       "depassement_pct": o["depassement_pct"], "e_max_V": o["e_max_V"],
                                       "fenetres_gains_changes": nb, "x_final": x_fin,
                                       "evenements": [{"t_ms": e_["t_ms"], "IAE_mVs": e_["IAE"] * 1e3,
                                                       "e_max_V": e_["e_max_V"], "revenu": e_["revenu"]} for e_ in ev]}
        print(f"    {nom:22s} apres 30 ms {o['IAE'] * 1e3:7.2f} ; demarrage {r['IAE_dem'] * 1e3:7.2f} ; "
              f"depassement {o['depassement_pct']:5.1f} % ; ecart max {o['e_max_V']:5.2f} V ; gains changes {nb:2d} fois")
        for e_ in ev:
            print(f"        evenement a {e_['t_ms']:5.1f} ms : IAE {e_['IAE'] * 1e3:6.2f} mV.s, ecart max {e_['e_max_V']:5.2f} V")

with open(os.path.join(DOSSIER, "loi_incrementale_contre_bloc_commun.json"), "w", encoding="utf-8") as f:
    json.dump(SORTIE, f, indent=1)

# 3. Figure : S4 (deux evenements), S2 (saut de consigne), S1 (demarrage)
COUL = {"Ziegler-Nichols": "0.55", "incrementale (Lu)": "tab:purple", "bloc PID commun": "tab:blue"}
FENETRES = (("S4", 49.0, 56.0, "S4 : 50 ms, 5 -> 4 ohms"), ("S4", 69.0, 76.0, "S4 : 70 ms, 4 -> 5 ohms"),
            ("S2", 49.0, 56.0, "S2 : 50 ms, saut de consigne F1"), ("S1", 0.0, 10.0, "S1 : demarrage"))
fig, ax = plt.subplots(3, 4, figsize=(20, 10), sharex="col", gridspec_kw={"height_ratios": [2, 1, 1]})
for col, (code, t0, t1, titre_) in enumerate(FENETRES):
    for nom, r in RES[code].items():
        s = r["sim"]
        t = s["t"] * 1e3
        k = (t >= t0) & (t <= t1)
        lab = f"{nom} (IAE {r['o']['IAE'] * 1e3:.2f} mV.s apres 30 ms)" if code != "S1" else \
            f"{nom} (IAE {r['IAE_dem'] * 1e3:.1f} mV.s de 0 a 30 ms)"
        ax[0, col].plot(t[k], s["v"][k], color=COUL[nom], lw=1.0, label=lab)
        ax[1, col].plot(t[k], s["d"][k], color=COUL[nom], lw=0.7)
        if hasattr(r["reg"], "K_pas"):
            ns = FORMES[nom]
            Kx = np.array(r["reg"].K_pas) / ns["K_ZN"]
            ax[2, col].plot(t[k], Kx[k, 0], color=COUL[nom], lw=1.2)
    ax[0, col].plot(t[k], RES[code]["Ziegler-Nichols"]["sim"]["consigne"][k], "k--", lw=0.6)
    ax[0, col].set_title(titre_, fontsize=10, loc="left")
    if code == "S1":
        ax[0, col].set_ylim(0, 125)
    ax[0, col].legend(fontsize=7, loc="lower right")
    ax[2, col].set_xlabel("temps (ms)")
    for a in ax[:, col]:
        a.grid(alpha=0.3)
ax[0, 0].set_ylabel("Vout (V)")
ax[1, 0].set_ylabel("rapport cyclique")
ax[2, 0].set_ylabel("gain proportionnel\n(multiple de ZN)")
fig.suptitle("ELM-PID : loi incrementale de Lu contre bloc PID commun, memes reglages d'adaptation",
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER, "loi_incrementale_contre_bloc_commun.png"), dpi=110)
plt.close(fig)
print("\nloi_incrementale_contre_bloc_commun.json et loi_incrementale_contre_bloc_commun.png ecrits.")
