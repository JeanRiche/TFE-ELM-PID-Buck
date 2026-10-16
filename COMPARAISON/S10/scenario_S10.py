# =============================================================================
# scenario_S10.py
#
# OBJECTIF
#   Simuler le scenario S10 (changement durable du point de fonctionnement :
#   Vin = 160 V et charge 25 ohms a partir de 50 ms, puis echelons de charge
#   de 15 ms au nouveau point) avec les cinq methodes, chacune prise dans
#   son dossier :
#     1. construction de S10 et controle qu'il n'est aucun des onze essais ;
#     2. cinq methodes au nominal, puis L ou C a +-0.1 % (robustesse) ;
#     3. ELM-PID et PINN-PID avec l'adaptation coupee a 50 ms (gains figes),
#        et controle : coupure a 200 ms = forme adaptative, a l'identique ;
#     4. tableau console, classement des quatre methodes avec egalites (5 %) ;
#     5. resultats_S10.json et S10_vout_commande_gains.png.
#   Aucune methode n'est modifiee : les classes sont executees depuis leurs
#   fichiers, la coupure passe seulement le drapeau existant adapt.adapter a
#   False.
#
# DOSSIERS NECESSAIRES
#   ../../PSO_PID, ../../FUZZY_PID, ../../ELM_PID, ../../PINN_PID
#
# COMMENT LANCER CE SCRIPT
#   python scenario_S10.py
#   (quatre processus ; quelques minutes)
# =============================================================================

import os
import json
import time
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
RACINE = os.path.dirname(os.path.dirname(DOSSIER))
T_DEBUT = time.time()


def charger(dossier, fichier, marque):
    """Execute les definitions d'un banc (jusqu'a la marque) dans un espace
    de noms propre a la methode ."""
    chemin = os.path.join(RACINE, dossier, fichier)
    ns = {"__file__": chemin, "__name__": "definitions"}
    with open(chemin, encoding="utf-8") as f:
        src = f.read()
    exec(src[:src.index(marque)], ns)
    return ns


NS = {"Ziegler-Nichols": charger("ELM_PID", "banc_elm_pid.py", "# %% ETAPE 2 :"),
      "PSO-PID": charger("PSO_PID", "recherche_pso_pid.py", "\n# %% ETAPE 3 :"),
      "Fuzzy-PID": charger("FUZZY_PID", "banc_fuzzy_pid.py", "# %% ETAPE 2 :"),
      "ELM-PID": None, "PINN-PID": charger("PINN_PID", "banc_pinn_pid.py", "# %% ETAPE 2 :")}
NS["ELM-PID"] = NS["Ziegler-Nichols"]
with open(os.path.join(RACINE, "PSO_PID", "predictions_banc_pso_pid.json"), encoding="utf-8") as f:
    G_PSO = json.load(f)["gains"]
FABRIQUES = {"Ziegler-Nichols": lambda ns: ns["PIDClassique"](),
             "PSO-PID": lambda ns: ns["PIDParallele"](G_PSO["P"], G_PSO["I"], G_PSO["D"]),
             "Fuzzy-PID": lambda ns: ns["FuzzyPID"](),
             "ELM-PID": lambda ns: ns["ELMPID"](),
             "PINN-PID": lambda ns: ns["PINNPID"]()}
METHODES = list(FABRIQUES)
CLASSEES = ["PSO-PID", "Fuzzy-PID", "ELM-PID", "PINN-PID"]
ADAPTATIVES = ["ELM-PID", "PINN-PID"]
NSB = NS["ELM-PID"]                           # constantes du banc commun (identiques dans les quatre dossiers)
TC = NSB["TC"]
K_ZN = np.array([NSB["PID_P"], NSB["PID_I"], NSB["PID_D"]])
for m in METHODES:
    if NS[m]["TC"] != TC or NS[m]["L_BOB"] != NSB["L_BOB"] or NS[m]["C_CONV"] != NSB["C_CONV"]:
        raise RuntimeError(f"{m} : banc commun different (Tc, L ou C).")


# %% ETAPE 1 : construction de S10 et controle d'identite

def construire_S10():
    """Construction de S10 (meme construction que les essais de reglage de PSO_PID/recherche_pso_pid.py)."""
    duree = 0.200
    n = int(round(duree / TC)) + 1
    t = np.arange(n) * TC

    def palier(base, morceaux):
        y = np.full(n, base, dtype=float)
        for t0, t1, val in morceaux:
            y[(t >= t0 - 1e-12) & (t < t1 - 1e-12)] = val
        return y

    R0 = 5.0
    g = lambda R: 1.0 / R - 1.0 / R0
    gx = palier(0.0, [(0.050, 1.0, g(25.0)), (0.100, 0.115, g(20.0)), (0.130, 0.145, g(30.0)),
                      (0.160, 0.175, g(20.0))])
    return {"code": "S10", "nom": "Changement durable du point", "duree": duree, "R0": R0, "q": 1e-9,
            "evenements": [0.05, 0.10, 0.115, 0.13, 0.145, 0.16, 0.175], "t": t,
            "dvref": np.zeros(n), "vin": palier(200.0, [(0.050, 1.0, 160.0)]), "gx": gx,
            "bruit": np.zeros(n)}


S10 = construire_S10()
N = len(S10["t"])
K50, K100, K200 = int(round(0.05 / TC)), int(round(0.10 / TC)), N - 1
print(f"S10 : {N} instants Tc, k(50 ms) = {K50}, k(100 ms) = {K100}, dernier k = {K200}.")
# Charge totale et Vin aux instants caracteristiques (controle de construction)
for tms in (10, 60, 105, 120, 135, 150, 165, 190):
    k = int(round(tms * 1e-3 / TC))
    print(f"  t = {tms:3d} ms : Vin = {S10['vin'][k]:5.1f} V, charge = {1.0 / (1.0 / S10['R0'] + S10['gx'][k]):5.2f} ohms")
for s in NSB["SCENARIOS"]:
    meme = (len(s["t"]) == N and abs(s["R0"] - S10["R0"]) < 1e-12 and np.array_equal(s["vin"], S10["vin"])
            and np.array_equal(s["gx"], S10["gx"]) and np.array_equal(s["dvref"], S10["dvref"]))
    if meme:
        raise RuntimeError(f"S10 coincide avec l'essai {s['code']}.")
print("  Controle : S10 ne coincide avec aucun des onze essais (R0, Vin, gx, dvref compares).")


# %% ETAPE 2 : simulations (un processus par simulation)

class Coupure:
    """Enrobage minimal : passe le drapeau existant adapt.adapter a False
    juste avant le pas k_coupure ; tout le reste est le regulateur lui-meme."""
    utilise_mesure = True

    def __init__(self, reg, k_coupure):
        self.reg, self.kc, self.k = reg, k_coupure, 0
        self.K_pas = reg.K_pas
        self.journal = reg.journal

    def pas(self, e, mesure):
        if self.k == self.kc:
            self.reg.adapt.adapter = False
        self.k += 1
        return self.reg.pas(e, mesure)


VARIANTES = (("nominal", 1.0, 1.0), ("L - 0.1 %", 0.999, 1.0), ("L + 0.1 %", 1.001, 1.0),
             ("C - 0.1 %", 1.0, 0.999), ("C + 0.1 %", 1.0, 1.001))
TACHES = [(m, nom, fl, fc, None) for m in METHODES for nom, fl, fc in VARIANTES]
TACHES += [(m, "fige a 50 ms", 1.0, 1.0, K50) for m in ADAPTATIVES]
TACHES += [(m, "coupure a 200 ms", 1.0, 1.0, K200) for m in ADAPTATIVES]


def une_simulation(tache):
    m, nom, fl, fc, kc = tache
    ns = NS[m]
    t0 = time.time()
    L0, C0 = ns["L_BOB"], ns["C_CONV"]
    ns["L_BOB"], ns["C_CONV"] = L0 * fl, C0 * fc             # le circuit seul change
    try:
        reg = FABRIQUES[m](ns)
        if kc is not None:
            if not hasattr(reg, "adapt") or not hasattr(reg.adapt, "adapter"):
                raise RuntimeError(f"{m} : pas de drapeau adapt.adapter.")
            reg = Coupure(reg, kc)
        sim = ns["simuler"](reg, S10)
        o = ns["grandeurs"](sim, S10)
    finally:
        ns["L_BOB"], ns["C_CONV"] = L0, C0
    e = sim["consigne"] - sim["v"]
    n = np.arange(N)
    w = n >= K100
    o["IAE_100_fin"] = float(np.sum(np.abs(e[w])) * TC)
    o["e_max_100_fin_V"] = float(np.max(np.abs(e[w])))
    if hasattr(reg, "K_pas") and len(reg.K_pas):
        K = np.array(reg.K_pas, dtype=float)
    else:
        P_ = getattr(reg, "P", NSB["PID_P"]); I_ = getattr(reg, "I", NSB["PID_I"]); D_ = getattr(reg, "D", NSB["PID_D"])
        K = np.tile([P_, I_, D_], (N, 1))
    sig = {"v": sim["v"], "d": sim["d"], "iL": sim["iL"], "K": K} if (fl, fc) == (1.0, 1.0) else None
    return (m, nom), {"o": o, "sig": sig, "duree_s": time.time() - t0}


if __name__ == "__main__":
    with mp.get_context("fork").Pool(4) as pool:
        RES = dict(pool.map(une_simulation, TACHES, chunksize=1))
    print(f"\n{len(TACHES)} simulations en {time.time() - T_DEBUT:.0f} s.")

    IAE = {m: {nom: RES[(m, nom)]["o"]["IAE_100_fin"] * 1e3 for nom, _, _ in VARIANTES} for m in METHODES}
    nomi = {m: RES[(m, "nominal")]["o"] for m in METHODES}

    # %% ETAPE 3 : controles de la coupure de l'adaptation
    titre = lambda txt: print("\n" + "=" * 100 + "\n " + txt + "\n" + "=" * 100)
    titre("ETAPE 3 : controles de la variante gains figes")
    for m in ADAPTATIVES:
        a, c200, c50 = RES[(m, "nominal")]["sig"], RES[(m, "coupure a 200 ms")]["sig"], RES[(m, "fige a 50 ms")]["sig"]
        identique = all(np.array_equal(a[x], c200[x]) for x in ("v", "d", "K"))
        avant = (np.array_equal(a["v"][:K50 + 1], c50["v"][:K50 + 1]) and np.array_equal(a["d"][:K50], c50["d"][:K50])
                 and np.array_equal(a["K"][:K50 + 1], c50["K"][:K50 + 1]))
        constants = bool(np.all(c50["K"][K50:] == c50["K"][K50]))
        print(f"  {m} : coupure a 200 ms = adaptatif a l'identique (Vout, d, gains) : {identique} ; "
              f"fige a 50 ms identique avant k = {K50} : {avant} ; gains constants de k = {K50} a la fin : {constants}")
        print(f"      multiplicateurs de ZN figes a 50 ms : ({', '.join(f'{x:.3f}' for x in c50['K'][K50] / K_ZN)}) ; "
              f"adaptatif a 100 ms : ({', '.join(f'{x:.3f}' for x in a['K'][K100] / K_ZN)}) ; "
              f"a la fin : ({', '.join(f'{x:.3f}' for x in a['K'][-1] / K_ZN)})")
        if not (identique and avant and constants):
            raise RuntimeError(f"{m} : la coupure de l'adaptation ne se comporte pas comme prevu.")

    # %% ETAPE 4 : tableau, classement, robustesse, adaptation
    titre("ETAPE 4a : IAE de 100 ms a la fin (mV.s), nominal ; grandeurs publiees a cote")
    print(f"  {'methode':18s} {'IAE 100-fin':>11s} {'/ ZN':>6s} {'emax 100-fin':>12s} {'IAE 50-100':>10s} "
          f"{'emax a 50 ms':>12s} {'DCM %':>6s}")
    lignes = METHODES + [f"{m} fige" for m in ADAPTATIVES]
    O = dict(nomi, **{f"{m} fige": RES[(m, "fige a 50 ms")]["o"] for m in ADAPTATIVES})
    zn = O["Ziegler-Nichols"]["IAE_100_fin"]
    for m in lignes:
        o = O[m]
        print(f"  {m:18s} {o['IAE_100_fin'] * 1e3:11.2f} {o['IAE_100_fin'] / zn:6.3f} {o['e_max_100_fin_V']:12.2f} "
              f"{o['evenements'][0]['IAE'] * 1e3:10.2f} {o['evenements'][0]['e_max_V']:12.2f} {o['dcm_pct']:6.2f}")

    titre("ETAPE 4b : par echelon (IAE mV.s / ecart max V / retour dans +-1 V)")
    NOMS_EV = ["100 ms 25->20", "115 ms 20->25", "130 ms 25->30", "145 ms 30->25", "160 ms 25->20", "175 ms 20->25"]
    print("  " + " " * 18 + "".join(f"{x:>24s}" for x in NOMS_EV))

    def retour(ev):
        if ev["reste_dans_bande"]:
            return "dans bande"
        return f"{ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU"

    for m in lignes:
        evs = O[m]["evenements"][1:]
        print(f"  {m:18s}" + "".join(f"{ev['IAE'] * 1e3:7.2f} /{ev['e_max_V']:5.2f} /{retour(ev):>10s}" for ev in evs))

    def classement(valeurs):
        o = sorted(valeurs, key=valeurs.get)
        txt = o[0]
        for a, b in zip(o, o[1:]):
            txt += (" = " if valeurs[b] / valeurs[a] - 1 < 0.05 else " < ") + b
        return o, txt

    titre("ETAPE 4c : classement des quatre methodes (IAE de 100 ms a la fin ; '=' : moins de 5 %)")
    ORDRES = {}
    for nom, _, _ in VARIANTES:
        val = {m: IAE[m][nom] for m in CLASSEES}
        o, txt = classement(val)
        ORDRES[nom] = o
        ecarts = ", ".join(f"{b}/{a} {100 * (val[b] / val[a] - 1):.1f} %" for a, b in zip(o, o[1:]))
        print(f"  {nom:10s} {txt}\n             ({ecarts})")
    print("\n  IAE par variante du circuit (mV.s) :")
    for m in METHODES:
        print(f"  {m:16s} " + "  ".join(f"{n} {IAE[m][n]:6.2f}" for n, _, _ in VARIANTES))
    paires = [(a, b) for i, a in enumerate(ORDRES["nominal"]) for b in ORDRES["nominal"][i + 1:]
              if IAE[b]["nominal"] / IAE[a]["nominal"] - 1 >= 0.05]
    retournements = [(nom, a, b) for nom, _, _ in VARIANTES for a, b in paires if IAE[a][nom] > IAE[b][nom]]
    robuste = not retournements
    print(f"  Robustesse : {len(paires)} paire(s) separee(s) d'au moins 5 % au nominal ; "
          + ("aucune ne s'inverse : classement robuste." if robuste else f"retournements : {retournements}"))

    titre("ETAPE 4d : contribution de l'adaptation (IAE de 100 ms a la fin, mV.s)")
    GAIN = {}
    for m in ADAPTATIVES:
        a, f_ = nomi[m]["IAE_100_fin"], RES[(m, "fige a 50 ms")]["o"]["IAE_100_fin"]
        GAIN[m] = 1.0 - a / f_
        verdict = ("sans effet mesurable (moins de 5 %)" if abs(f_ / a - 1) < 0.05 else
                   ("l'adaptation aide" if a < f_ else "l'adaptation NUIT"))
        print(f"  {m} : adaptatif {a * 1e3:.2f}, fige a 50 ms {f_ * 1e3:.2f} ; gain 1 - adaptatif/fige = "
              f"{100 * GAIN[m]:+.1f} % : {verdict}")
        print(f"      IAE 50-100 ms : adaptatif {nomi[m]['evenements'][0]['IAE'] * 1e3:.2f}, fige "
              f"{RES[(m, 'fige a 50 ms')]['o']['evenements'][0]['IAE'] * 1e3:.2f}")

    v = {m: IAE[m]["nominal"] for m in METHODES}

    # %% ETAPE 5 : fichiers
    resultats = {"scenario": "S10",
                 "definition": {"R0": S10["R0"], "evenements": S10["evenements"], "duree": S10["duree"],
                                "q": S10["q"], "Vin_apres_50ms": 160.0,
                                "charges": "25 ohms ; 20 sur [100;115), 30 sur [130;145), 20 sur [160;175) ms"},
                 "IAE_100_fin_mVs": IAE,
                 "IAE_100_fin_fige_50ms_mVs": {m: RES[(m, "fige a 50 ms")]["o"]["IAE_100_fin"] * 1e3 for m in ADAPTATIVES},
                 "gain_adaptation": GAIN,
                 "ordres": ORDRES, "classement_nominal": classement({m: v[m] for m in CLASSEES})[1],
                 "robuste": robuste, "retournements": retournements,
                 "grandeurs_nominal": {m: O[m] for m in lignes},
                 "multiplicateurs_ZN": {m: {"50ms": (RES[(m, "nominal")]["sig"]["K"][K50] / K_ZN).tolist(),
                                            "100ms": (RES[(m, "nominal")]["sig"]["K"][K100] / K_ZN).tolist(),
                                            "fin": (RES[(m, "nominal")]["sig"]["K"][-1] / K_ZN).tolist()}
                                        for m in ADAPTATIVES + ["Fuzzy-PID"]},
                 "duree_s": time.time() - T_DEBUT}
    with open(os.path.join(DOSSIER, "resultats_S10.json"), "w", encoding="utf-8") as f:
        json.dump(resultats, f, indent=1)

    COUL = {"Ziegler-Nichols": "0.5", "PSO-PID": "tab:orange", "Fuzzy-PID": "tab:green", "ELM-PID": "tab:blue",
            "PINN-PID": "tab:red"}
    t_ms = S10["t"] * 1e3
    k = t_ms >= 45.0
    fig, ax = plt.subplots(5, 1, figsize=(14, 15), sharex=True, gridspec_kw={"height_ratios": [2.2, 1, 1, 1, 1]})
    courbes = [(m, RES[(m, "nominal")]["sig"], "-", f"{m} (IAE 100-200 ms {v[m]:.2f} mV.s)") for m in METHODES]
    courbes += [(m, RES[(m, "fige a 50 ms")]["sig"], "--",
                 f"{m} fige a 50 ms ({RES[(m, 'fige a 50 ms')]['o']['IAE_100_fin'] * 1e3:.2f} mV.s)") for m in ADAPTATIVES]
    for m, s, st, lab in courbes:
        ax[0].plot(t_ms[k], s["v"][k], st, color=COUL[m], lw=0.8, label=lab)
        ax[1].plot(t_ms[k], s["d"][k], st, color=COUL[m], lw=0.6)
        for j in range(3):
            ax[2 + j].plot(t_ms[k], s["K"][k, j] / K_ZN[j], st, color=COUL[m], lw=1.0)
    ax[0].axhline(100.0, color="k", lw=0.5)
    ax[0].axhspan(99.0, 101.0, color="0.9", zorder=0)
    ax[0].set_ylim(96.0, 104.0)
    ax[0].set_ylabel("Vout (V)\n(zoom, bande +-1 V grisee)")
    ax[1].set_ylabel("rapport cyclique")
    ax[1].set_ylim(0.45, 0.75)
    for j, g in enumerate(("P", "I", "D")):
        ax[2 + j].set_ylabel(f"{g} / {g} de ZN")
    for a in ax:
        a.grid(alpha=0.3)
        for te in S10["evenements"]:
            a.axvline(te * 1e3, color="k", lw=0.4, ls=":")
    ax[0].axvspan(50.0, 100.0, color="0.97", zorder=0)
    ax[0].text(51.0, 103.4, "50 ms : Vin 200 -> 160 V, charge 5 -> 25 ohms ; reajustement non mesure", fontsize=8)
    for te, txt in ((100, "20 ohms"), (130, "30 ohms"), (160, "20 ohms")):
        ax[0].text(te + 0.5, 96.4, txt, fontsize=8)
    ax[0].legend(fontsize=8, loc="upper right", ncol=2)
    ax[-1].set_xlabel("temps (ms)")
    ax[-1].set_xlim(45.0, 200.0)
    fig.suptitle("S10 : changement durable du point de fonctionnement, puis echelons de charge au nouveau point "
                 "(banc commun)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(DOSSIER, "S10_vout_commande_gains.png"), dpi=110)
    plt.close(fig)
    print(f"\n  resultats_S10.json et S10_vout_commande_gains.png ecrits. Duree totale {time.time() - T_DEBUT:.0f} s.")
