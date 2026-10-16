# =============================================================================
# banc_pso_pid.py
#
# OBJECTIF
#   Juger sur les onze essais de developpement le reglage trouve par
#   recherche_pso_pid.py, a cote de Ziegler-Nichols
#   (meme PID, gains de depart) et du meilleur PID fige sur grille, et ecrire
#   ce que Simulink doit retrouver :
#     - predictions_banc_pso_pid.json : gains retenus et resultats attendus
#       (meme format que predictions_banc_pid_fige.json, lu par
#       Construction_PSO_PID.m et Simuler_PSO_PID.m) ;
#     - banc_pso_pid_resultats.json : tous les chiffres ;
#     - banc_pso_pid.png (S1, S5, S8b, S9) et pso_pid_convergence.png.
#   Les onze essais n'ont servi a rien pendant la recherche.
#   Retenu : la graine de plus petit J_reglage.
#
# COMMENT LANCER CE SCRIPT
#   python banc_pso_pid.py      (une minute)
# =============================================================================


# %% ETAPE 0 : definitions et recherches

import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt

try:
    DOSSIER_PSO = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PSO = os.getcwd()
with open(os.path.join(DOSSIER_PSO, "recherche_pso_pid.py"), encoding="utf-8") as f:
    _SRC = f.read()
exec(_SRC[:_SRC.index("\n# %% ETAPE 3 :")])   # banc commun, contraintes, essais E1 a E4, J_reglage

T_DEBUT = time.time()
ESSAIS = {sc["code"]: sc for sc in SCENARIOS}
CODES_IAE = ["S1", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]
CODES_DEM = ["S1", "S8a", "S8b"]
K30_ESSAI = int(round(0.03 / TC))
with open(FICHIER, encoding="utf-8") as f:
    RECH = json.load(f)["recherches"]
CLES = [f"graine={g}" for g in GRAINES]
manque = [c for c in CLES if c not in RECH]
if manque:
    raise RuntimeError("Recherches manquantes dans recherche_pso_pid.json : " + ", ".join(manque))
with open(os.path.join(DOSSIER_PSO, "predictions_banc_pid_fige.json"), encoding="utf-8") as f:
    GAINS_FIGE = json.load(f)["gains"]


def titre(texte):
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78)


# %% ETAPE 1 : les cinq recherches et le reglage retenu

titre("ETAPE 1 : les cinq recherches et le reglage retenu")
print(f"  {'recherche':10s} {'J_reglage':>10s} {'a (P)':>8s} {'b (I)':>8s} {'c (D)':>8s} {'marge':>7s} {'coupure':>8s}")
for c in CLES:
    r = RECH[c]
    print(f"  {c:10s} {r['J_reglage']:10.5f} {r['gbest'][0]:8.4f} {r['gbest'][1]:8.4f} {r['gbest'][2]:8.4f} "
          f"{r['marge_min']:6.1f}d {r['fc_max']:7.0f}Hz")
CLE_RETENUE = min((c for c in CLES if RECH[c]["admissible"]), key=lambda c: RECH[c]["J_reglage"])
RET = RECH[CLE_RETENUE]
P_R, I_R, D_R = RET["gains"]["P"], RET["gains"]["I"], RET["gains"]["D"]
print(f"\n  Retenu : {CLE_RETENUE} ; P = {P_R:.6f}, I = {I_R:.4f}, D = {D_R:.6e} "
      f"(x ZN : {RET['gbest'][0]:.4f}, {RET['gbest'][1]:.4f}, {RET['gbest'][2]:.4f}).")
print("  Rapports a Ziegler-Nichols sur les essais de reglage (demarrage / apres 30 ms) :")
for code, q in RET["ratios"].items():
    print(f"    {code} : {q['demarrage']:.3f} / {q['apres_30ms']:.3f}")
for nom, g in (("meilleur PID fige", GAINS_FIGE),):
    x = (g["P"] / PID_P, g["I"] / PID_I, g["D"] / PID_D)
    r = evaluer(x)
    print(f"  Pour information, J_reglage du {nom} : " +
          (f"{r['J']:.4f}" if r["admissible"] else f"refuse (violation {r['violation']:.3f})"))


# %% ETAPE 2 : jugement sur les onze essais

def evaluer_essais(P, I, D, essais=None):
    res = {}
    for code, sc in (essais or ESSAIS).items():
        sim = simuler(PIDParallele(P, I, D), sc)
        o = grandeurs(sim, sc)
        e = sim["consigne"] - sim["v"]
        res[code] = {"o": o, "IAE_dem": float(np.sum(np.abs(e[:K30_ESSAI])) * TC),
                     "revenus": bool(all(ev["revenu"] for ev in o["evenements"])), "sim": sim}
    return res


def termes(res):
    return np.array([res[c]["o"]["IAE"] for c in CODES_IAE] + [res[c]["IAE_dem"] for c in CODES_DEM])


titre("ETAPE 2 : jugement sur les onze essais communs")
TOUS = {"Ziegler-Nichols": evaluer_essais(PID_P, PID_I, PID_D),
        "meilleur PID fige": evaluer_essais(GAINS_FIGE["P"], GAINS_FIGE["I"], GAINS_FIGE["D"]),
        "PSO-PID": evaluer_essais(P_R, I_R, D_R)}
T_ZN = termes(TOUS["Ziegler-Nichols"])
J_TOUS = {n: float(np.mean(termes(r) / T_ZN)) for n, r in TOUS.items()}
print("  " + " " * 28 + "".join(f"{n:>19s}" for n in TOUS))
print("  " + f"{'cout J':28s}" + "".join(f"{J_TOUS[n]:19.3f}" for n in TOUS))
print("  " + f"{'evenements non revenus':28s}" + "".join(
    f"{(','.join(c for c in CODES_IAE if not TOUS[n][c]['revenus']) or '-'):>19s}" for n in TOUS))
for code in ESSAIS:
    for cle, etiquette, fmt, k in (("depassement_pct", "depassement %", "{:.1f}", 1.0), ("IAE", "IAE mV.s", "{:.2f}", 1e3),
                                   ("e_eff_V", "ecart efficace mV", "{:.1f}", 1e3), ("iL_max_dem_A", "iL crete dem. A", "{:.1f}", 1.0)):
        if cle == "depassement_pct" and code not in ("S1", "S8a", "S8b"):
            continue
        if cle == "iL_max_dem_A" and code not in ("S1", "S8b"):
            continue
        print("  " + f"{code + ' ' + etiquette:28s}" + "".join(f"{fmt.format(TOUS[n][code]['o'][cle] * k):>19s}" for n in TOUS))
print("  " + f"{'activite |dd| (S1)':28s}" + "".join(f"{TOUS[n]['S1']['o']['dd_moyen']:19.4f}" for n in TOUS))
DECOMP = {}
print("\n  Decomposition de J (10 IAE apres 30 ms ; 3 IAE de demarrage) :")
for n in TOUS:
    r = termes(TOUS[n]) / T_ZN
    DECOMP[n] = (float(np.mean(r[:10])), float(np.mean(r[10:])))
    print(f"    {n:20s} apres 30 ms {DECOMP[n][0]:.3f} ; demarrage {DECOMP[n][1]:.3f}")
print("\n  Rapports a ZN, terme par terme (" + ", ".join(CODES_IAE) + " ; demarrages " + ", ".join(CODES_DEM) + ") :")
for n in TOUS:
    print(f"    {n:20s} " + " ".join(f"{x:5.2f}" for x in termes(TOUS[n]) / T_ZN))

print("\n  Les cinq recherches jugees sur les onze essais (retenue : la graine de plus petit J_reglage) :")
JUGEMENT = {}
for c in CLES:
    r = RECH[c]
    res = TOUS["PSO-PID"] if c == CLE_RETENUE else evaluer_essais(r["gains"]["P"], r["gains"]["I"], r["gains"]["D"])
    Jc = float(np.mean(termes(res) / T_ZN))
    non = [k for k in CODES_IAE if not res[k]["revenus"]]
    JUGEMENT[c] = {"J": Jc, "non_revenus": non}
    print(f"    {c:10s} J = {Jc:.3f} ; {'tous revenus' if not non else 'pas revenus : ' + ','.join(non)}"
          f"{'   <- retenu' if c == CLE_RETENUE else ''}")


# %% ETAPE 3 : reglage retenu

PSO = TOUS["PSO-PID"]
J_pso = J_TOUS["PSO-PID"]


# %% ETAPE 3b : essais complementaires, hors du cout J
# S10 est simule ici avec les memes gains, pour figurer comme les onze essais dans
# predictions_banc_pso_pid.json (lu par Simuler_PSO_PID.m) et dans
# banc_pso_pid_resultats.json. Il n'entre pas dans J.

ESSAIS_HORS_J = {sc["code"]: sc for sc in SCENARIOS_COMPLEMENTAIRES}
HORS_J = {}
if ESSAIS_HORS_J:
    titre("ETAPE 3b : essais complementaires, hors du cout J (" + ", ".join(ESSAIS_HORS_J) + ")")
    HORS_J = {"Ziegler-Nichols": evaluer_essais(PID_P, PID_I, PID_D, ESSAIS_HORS_J),
              "PSO-PID": evaluer_essais(P_R, I_R, D_R, ESSAIS_HORS_J)}
    for code in ESSAIS_HORS_J:
        for n, res in HORS_J.items():
            print(f"  {code} {n:16s} {resume(res[code]['o'])}")
            for ev in res[code]["o"]["evenements"]:
                etat = ("reste dans la bande" if ev["reste_dans_bande"] else
                        (f"retour en {ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU"))
                print(f"      evenement a {ev['t_ms']:6.1f} ms : IAE {ev['IAE'] * 1e3:7.2f} mV.s, ecart max "
                      f"{ev['e_max_V']:6.2f} V, {etat}")


# %% ETAPE 4 : fichiers pour Simulink, resultats et figures

titre("ETAPE 4 : fichiers pour Simulink, resultats et figures")
pas_1ms = int(round(1e-3 / TC))
predictions = {"format": 2, "Te": TC,
               "regulateur": f"PSO-PID ({CLE_RETENUE} ; P = {P_R:.6f}, I = {I_R:.4f}, D = {D_R:.6e}, N = {PID_N})",
               "gains": {"P": P_R, "I": I_R, "D": D_R, "N": PID_N, "multiplicateurs": RET["gbest"]},
               "critere": {"J_reglage": RET["J_reglage"], "J_essais": J_pso, "marge_min_deg": RET["marge_min"],
                           "fc_max_Hz": RET["fc_max"]},
               "controle": {}, "essais": {}}
for code, r in list(PSO.items()) + list(HORS_J.get("PSO-PID", {}).items()):   # onze essais, puis S10
    sim = r["sim"]
    predictions["essais"][code] = {"grandeurs": r["o"],
                                   "v_toutes_les_ms": sim["v"][::pas_1ms].tolist(),
                                   "iL_toutes_les_ms": sim["iL"][::pas_1ms].tolist(),
                                   "d_moyen_par_ms": [float(np.mean(sim["d"][i:i + pas_1ms]))
                                                      for i in range(0, len(sim["d"]) - pas_1ms + 1, pas_1ms)]}
with open(os.path.join(DOSSIER_PSO, "predictions_banc_pso_pid.json"), "w", encoding="utf-8") as f:
    json.dump(predictions, f, indent=1, allow_nan=False)
print("  predictions_banc_pso_pid.json ecrit.")

sortie = {"retenu": CLE_RETENUE, "gains": predictions["gains"],
          "J": J_TOUS, "decomposition": DECOMP,
          "jugement_cinq_recherches": JUGEMENT,
          "references": {n: {c: {"grandeurs": TOUS[n][c]["o"], "IAE_dem": TOUS[n][c]["IAE_dem"],
                                 "revenus": TOUS[n][c]["revenus"]} for c in ESSAIS} for n in TOUS},
          "essais_hors_J": {"codes": list(ESSAIS_HORS_J),
                            "note": "essais complementaires (S10), hors du cout J",
                            "references": {n: {c: {"grandeurs": HORS_J[n][c]["o"], "IAE_dem": HORS_J[n][c]["IAE_dem"],
                                                   "revenus": HORS_J[n][c]["revenus"]} for c in ESSAIS_HORS_J}
                                           for n in HORS_J}},
          "duree_s": round(time.time() - T_DEBUT, 1)}
with open(os.path.join(DOSSIER_PSO, "banc_pso_pid_resultats.json"), "w", encoding="utf-8") as f:
    json.dump(sortie, f, indent=1, allow_nan=False, default=lambda x: None)
print("  banc_pso_pid_resultats.json ecrit.")

COULEURS = (("Ziegler-Nichols", "0.55"), ("meilleur PID fige", "tab:green"),
            ("PSO-PID", "tab:red"))
fig, axes = plt.subplots(2, 4, figsize=(20, 7.5), sharex="col")
for col, (code, xlim) in enumerate((("S1", (0, 15)), ("S5", (40, 140)), ("S8b", None), ("S9", None))):
    for nom, couleur in COULEURS:
        sim = TOUS[nom][code]["sim"]
        axes[0, col].plot(sim["t"] * 1e3, sim["v"], lw=0.6, color=couleur, label=nom)
        axes[1, col].plot(sim["t"] * 1e3, sim["iL"], lw=0.4, color=couleur)
    axes[0, col].plot(sim["t"] * 1e3, sim["consigne"], "k--", lw=0.6)
    axes[0, col].set_ylim((90, 110) if code == "S5" else (30, 200))
    axes[0, col].set_title(f"{code} : tension de sortie", fontsize=9, loc="left")
    axes[1, col].set_title(f"{code} : courant de la bobine (observe)", fontsize=9, loc="left")
    axes[1, col].set_xlabel("temps (ms)")
    if xlim:
        axes[1, col].set_xlim(*xlim)
    for ax in axes[:, col]:
        ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7)
axes[0, 0].set_ylabel("V")
axes[1, 0].set_ylabel("A")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_PSO, "banc_pso_pid.png"), dpi=110)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for c in CLES:
    h = RECH[c]["historique"]
    axes[0].plot([x["iter"] for x in h], [x["J_gbest"] if x["J_gbest"] is not None else np.nan for x in h], lw=1.0, label=c)
    axes[1].plot([x["iter"] for x in h], [x["part_admissible"] for x in h], lw=0.8)
axes[0].set_title("J_reglage du gbest", fontsize=9, loc="left")
axes[1].set_title("part d'individus admissibles", fontsize=9, loc="left")
for ax in axes:
    ax.set_xlabel("iteration")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_PSO, "pso_pid_convergence.png"), dpi=110)
plt.close(fig)
print(f"  banc_pso_pid.png et pso_pid_convergence.png ecrits. Duree totale : {time.time() - T_DEBUT:.0f} s.")
