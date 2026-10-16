# =============================================================================
# ablation_reseau_elm.py
#
# OBJECTIF
#   Mesurer le role du reseau ELM sur S2, S3, S8a et S10 :
#     1. trois formes : ELM-PID (reglages finaux), "Jacobien constant" et
#        "Sans OS-ELM", exactement les ablations de
#        ../../ELM_PID/banc_elm_pid.py (dictionnaire ABLATIONS, etape 4 :
#        memes arguments passes a ELMPID) ;
#     2. cinq variantes du circuit (nominal, L ou C a +-0.1 %), regulateurs
#        inchanges (seuls L_BOB et C_CONV changent) ;
#     3. controles prealables (sinon arret) : ELM-PID nominal =
#        metriques_banc.csv a 1e-9 relatif pres sur S2, S3, S8a, S10 ;
#        ablations nominales sur S2, S3, S8a = banc_elm_pid_resultats.json
#        (cle "ablations") a 1e-9 relatif pres ;
#     4. pour chaque paire (ELM-PID contre une ablation, par scenario) :
#        separee si le signe de la difference d'IAE est le meme dans les cinq
#        variantes ; ecart nominal, minimal et maximal en % ;
#     5. sur S10, multiplicateurs de Ziegler-Nichols (P, I, D) en fin d'essai
#        de chaque forme, au nominal.
#   IAE de classement (meme calcul que metriques_banc.py) : S2, S3, S8a de
#   30 ms a la fin ; S10 de 100 ms a la fin. S10 est construit comme dans
#   ../S10/scenario_S10.py (construire_S10).
#
# DOSSIERS NECESSAIRES
#   ../../ELM_PID (complet), ../metriques/metriques_banc.csv
#
# FICHIERS ECRITS (dans ce dossier)
#   ablation_resultats.json
#
# COMMENT LANCER CE SCRIPT
#   python ablation_reseau_elm.py > ablation_sortie_console.txt
#   (quatre processus ; plusieurs minutes)
# =============================================================================

import os
import json
import time
import multiprocessing as mp
import numpy as np

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
COMP = os.path.dirname(DOSSIER)
RACINE = os.path.dirname(COMP)
T_DEBUT = time.time()


def charger(dossier, fichier, marque):
    """Execute les definitions d'un banc (jusqu'a la marque) dans un espace
    de noms propre ."""
    chemin = os.path.join(RACINE, dossier, fichier)
    ns = {"__file__": chemin, "__name__": "definitions"}
    with open(chemin, encoding="utf-8") as f:
        src = f.read()
    exec(src[:src.index(marque)], ns)
    return ns, src


NS, SRC_ELM = charger("ELM_PID", "banc_elm_pid.py", "# %% ETAPE 2 :")

# Ablations : memes arguments que le dictionnaire ABLATIONS de banc_elm_pid.py
# (etape 4). Controle textuel : les deux lignes doivent figurer telles quelles
# dans le banc, sinon arret.
LIGNES_ABL = ['ABLATIONS = {"Jacobien constant": {"JACOBIEN_CONSTANT": J_REGRESSION},',
              '"Sans OS-ELM": {"APPRENDRE": 0.0},']
for l_ in LIGNES_ABL:
    if l_ not in SRC_ELM:
        raise RuntimeError(f"banc_elm_pid.py ne contient plus la definition attendue : {l_}")
FORMES = {"ELM-PID": {},
            "Jacobien constant": {"JACOBIEN_CONSTANT": NS["J_REGRESSION"]},
            "Sans OS-ELM": {"APPRENDRE": 0.0}}
NOMS_FORMES = list(FORMES)
ABL = ["Jacobien constant", "Sans OS-ELM"]
SCENARIOS = ["S2", "S3", "S8a", "S10"]
VARIANTES = (("nominal", 1.0, 1.0), ("L - 0.1 %", 0.999, 1.0), ("L + 0.1 %", 1.001, 1.0),
             ("C - 0.1 %", 1.0, 0.999), ("C + 0.1 %", 1.0, 1.001))
NOMS_V = [v[0] for v in VARIANTES]
TC = NS["TC"]
K_ZN = np.array([NS["PID_P"], NS["PID_I"], NS["PID_D"]])
K30 = int(round(0.03 / TC))
K100 = int(round(0.10 / TC))


def construire_S10():
    """Copie de construire_S10 de ../S10/scenario_S10.py (meme construction)."""
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


def scenario(code):
    if code == "S10":
        return S10
    return {s["code"]: s for s in NS["SCENARIOS"]}[code]


def fenetre_classement(code, N):
    """Meme fenetre que metriques_banc.py."""
    return {"S10": (K100, N)}.get(code, (K30, N))


def lancer(tache):
    """Comme separation_methodes.py : le circuit seul change (L_BOB, C_CONV)."""
    code, vers, nom, fl, fc = tache
    t0 = time.time()
    L0, C0 = NS["L_BOB"], NS["C_CONV"]
    NS["L_BOB"], NS["C_CONV"] = L0 * fl, C0 * fc
    try:
        sc = scenario(code)
        reg = NS["ELMPID"](**FORMES[vers])
        sim = NS["simuler"](reg, sc)
    finally:
        NS["L_BOB"], NS["C_CONV"] = L0, C0
    ka, kb = fenetre_classement(code, len(sim["v"]))
    iae = float(np.sum(np.abs((sim["consigne"] - sim["v"])[ka:kb])) * TC) * 1e3      # mV.s
    x_fin = (np.array(reg.K_pas[-1], dtype=float) / K_ZN).tolist()
    return (code, vers, nom), (iae, x_fin, time.time() - t0)


def lire_csv():
    ref = {}
    with open(os.path.join(COMP, "metriques", "metriques_banc.csv"), encoding="utf-8") as f:
        next(f)
        for ligne in f:
            s, m, g, fe, me, v, u = ligne.strip().split(";")
            if g == "classement" and fe == "classement" and me == "IAE":
                ref[(s, m)] = float(v)
    return ref


def signe(x):
    return (x > 0) - (x < 0)


if __name__ == "__main__":
    for s in NS.get("SCENARIOS_COMPLEMENTAIRES", []):
        if s["code"] == "S10":
            meme = all(np.array_equal(s[k], S10[k]) for k in ("t", "vin", "gx", "dvref")) and s["R0"] == S10["R0"]
            print(f"S10 construit comme scenario_S10.py ; identique a l'essai complementaire du banc : {meme}")
    print(f"Jacobien constant des ablations : {NS['J_REGRESSION']:.6f}")

    TACHES = [(c, v, nom, fl, fc) for c in SCENARIOS for v in NOMS_FORMES for nom, fl, fc in VARIANTES]
    TACHES.sort(key=lambda t: t[0] != "S10")
    print(f"{len(TACHES)} simulations, quatre processus.", flush=True)
    RES, XFIN = {}, {}
    with mp.get_context("fork").Pool(4) as pool:
        for i, (cle, (iae, x_fin, d)) in enumerate(pool.imap_unordered(lancer, TACHES, chunksize=1), 1):
            RES[cle], XFIN[cle] = iae, x_fin
            print(f"  [{i:3d}/{len(TACHES)}] {cle[0]:4s} {cle[1]:18s} {cle[2]:10s} IAE {iae:10.5f} mV.s "
                  f"({d:5.1f} s)", flush=True)
    print(f"Simulations terminees en {time.time() - T_DEBUT:.0f} s.")
    IAE = {c: {v: {nom: RES[(c, v, nom)] for nom in NOMS_V} for v in NOMS_FORMES} for c in SCENARIOS}

    # %% Controles prealables (tolerance 1e-9 relatif)
    REF = lire_csv()
    with open(os.path.join(RACINE, "ELM_PID", "banc_elm_pid_resultats.json"), encoding="utf-8") as f:
        REF_ABL = json.load(f)["ablations"]
    print("\nControles prealables (tolerance 1e-9 relatif)")
    ecart_max, echecs, controles = 0.0, [], []
    for c in SCENARIOS:
        controles.append((f"{c} ELM-PID / metriques_banc.csv", IAE[c]["ELM-PID"]["nominal"], REF[(c, "ELM-PID")]))
    for a in ABL:
        for c in ("S2", "S3", "S8a"):
            controles.append((f"{c} {a} / banc_elm_pid_resultats.json", IAE[c][a]["nominal"],
                              REF_ABL[a][c]["grandeurs"]["IAE"] * 1e3))
    for lib, x, r in controles:
        e = abs(x - r) / abs(r)
        ecart_max = max(ecart_max, e)
        print(f"  {lib:48s} {x:.10g} contre {r:.10g} (ecart {e:.2e})")
        if e > 1e-9:
            echecs.append(f"{lib} : {x:.10g} contre {r:.10g} (ecart {e:.2e})")
    print(f"  ecart relatif maximal {ecart_max:.2e}")
    if echecs:
        print("  ARRET : le calcul n'est pas valide.\n  " + "\n  ".join(echecs))
        with open(os.path.join(DOSSIER, "ablation_resultats.json"), "w", encoding="utf-8") as f:
            json.dump({"valide": False, "echecs_controle": echecs, "IAE_mVs": IAE}, f, indent=1)
        raise SystemExit(1)

    # %% Paires ELM-PID / ablation
    SORTIE = {"regle": "signe de la difference d'IAE identique dans les cinq variantes", "valide": True, "ecart_controle_max": ecart_max,
              "formes": {v: {k: float(x) for k, x in FORMES[v].items()} for v in NOMS_FORMES},
              "variantes": NOMS_V, "IAE_mVs": IAE, "scenarios": {}}
    for c in SCENARIOS:
        print("\n" + "=" * 100 + f"\n {c} : IAE de classement (mV.s) par variante\n" + "=" * 100)
        print("  " + " " * 18 + "".join(f"{n:>12s}" for n in NOMS_V))
        for v in NOMS_FORMES:
            print(f"  {v:18s}" + "".join(f"{IAE[c][v][n]:12.4f}" for n in NOMS_V))
        paires = []
        print(f"\n  Paires (ecart = IAE_ablation / IAE_ELM-PID - 1, en %)")
        for a in ABL:
            ec = {n: 100.0 * (IAE[c][a][n] / IAE[c]["ELM-PID"][n] - 1.0) for n in NOMS_V}
            sg = [signe(IAE[c][a][n] - IAE[c]["ELM-PID"][n]) for n in NOMS_V]
            sep = len(set(sg)) == 1 and sg[0] != 0
            paires.append({"A": "ELM-PID", "B": a, "separee": sep, "signes": sg,
                           "ecart_nominal_pct": ec["nominal"], "ecart_min_pct": min(ec.values()),
                           "ecart_max_pct": max(ec.values()), "ecart_pct": ec})
            print(f"    ELM-PID / {a:18s} {'SEPAREE    ' if sep else 'NON SEPAREE'}  nominal {ec['nominal']:+8.3f} %"
                  f"   min {min(ec.values()):+8.3f} %   max {max(ec.values()):+8.3f} %   signes {sg}")
        SORTIE["scenarios"][c] = {"paires": paires}

    SORTIE["multiplicateurs_ZN_fin_S10_nominal"] = {v: XFIN[("S10", v, "nominal")] for v in NOMS_FORMES}
    print("\n  S10, nominal : multiplicateurs de Ziegler-Nichols (P, I, D) en fin d'essai")
    for v in NOMS_FORMES:
        print(f"    {v:18s} (" + ", ".join(f"{x:.3f}" for x in XFIN[("S10", v, "nominal")]) + ")")

    SORTIE["duree_s"] = time.time() - T_DEBUT
    with open(os.path.join(DOSSIER, "ablation_resultats.json"), "w", encoding="utf-8") as f:
        json.dump(SORTIE, f, indent=1)
    print(f"\nablation_resultats.json ecrit. Duree totale {time.time() - T_DEBUT:.0f} s.")
