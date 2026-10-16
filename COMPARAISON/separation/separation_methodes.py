# =============================================================================
# separation_methodes.py
#
# OBJECTIF
#   Mesurer, sur chacun des cinq scenarios du corps du chapitre 5 (S1, S2,
#   S3, S8a, S10), quelles paires de methodes sont separees :
#     1. cinq methodes x cinq variantes du circuit (nominal, L ou C a
#        +-0.1 %), regulateurs inchanges (seuls L_BOB et
#        C_CONV changent) ;
#     2. controle : IAE nominale = metriques_banc.csv a 1e-9 relatif pres
#        (sinon arret) ;
#     3. pour chaque paire : separee si le signe de la difference d'IAE est
#        le meme dans les cinq variantes ; ecart nominal, minimal et maximal
#        en % ;
#     4. groupes ordonnes si la relation "non separees" est transitive,
#        sinon matrice des paires telle quelle.
#   IAE de classement (meme calcul que metriques_banc.py) : S1 de 0 a 30 ms ;
#   S2, S3, S8a de 30 ms a la fin ; S10 de 100 ms a la fin.
#   S10 est construit comme dans ../S10/scenario_S10.py (construire_S10).
#
# DOSSIERS NECESSAIRES
#   ../../PSO_PID, ../../FUZZY_PID, ../../ELM_PID, ../../PINN_PID (complets),
#   ../metriques/metriques_banc.csv
#
# FICHIERS ECRITS (dans ce dossier)
#   separation_resultats.json
#
# COMMENT LANCER CE SCRIPT
#   python separation_methodes.py > separation_sortie_console.txt
#   (quatre processus ; plusieurs minutes)
# =============================================================================

import os
import json
import time
import itertools
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
SCENARIOS = ["S1", "S2", "S3", "S8a", "S10"]
VARIANTES = (("nominal", 1.0, 1.0), ("L - 0.1 %", 0.999, 1.0), ("L + 0.1 %", 1.001, 1.0),
             ("C - 0.1 %", 1.0, 0.999), ("C + 0.1 %", 1.0, 1.001))
NOMS_V = [v[0] for v in VARIANTES]
NSB = NS["ELM-PID"]
TC = NSB["TC"]
for m_ in METHODES:
    if NS[m_]["TC"] != TC or NS[m_]["L_BOB"] != NSB["L_BOB"] or NS[m_]["C_CONV"] != NSB["C_CONV"]:
        raise RuntimeError(f"{m_} : banc commun different (Tc, L ou C).")
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


def scenario(ns, code):
    if code == "S10":
        return S10
    return {s["code"]: s for s in ns["SCENARIOS"]}[code]


def fenetre_classement(code, N):
    """Meme fenetre que metriques_banc.py."""
    return {"S1": (0, K30), "S10": (K100, N)}.get(code, (K30, N))


def lancer(tache):
    """Le circuit seul change (L_BOB, C_CONV)."""
    code, m, nom, fl, fc = tache
    ns = NS[m]
    t0 = time.time()
    L0, C0 = ns["L_BOB"], ns["C_CONV"]
    ns["L_BOB"], ns["C_CONV"] = L0 * fl, C0 * fc
    try:
        sc = scenario(ns, code)
        reg = FABRIQUES[m](ns)
        sim = ns["simuler"](reg, sc)
    finally:
        ns["L_BOB"], ns["C_CONV"] = L0, C0
    ka, kb = fenetre_classement(code, len(sim["v"]))
    iae = float(np.sum(np.abs((sim["consigne"] - sim["v"])[ka:kb])) * TC) * 1e3      # mV.s
    return (code, m, nom), (iae, time.time() - t0)


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
    # Controle de construction de S10 : identique a l'essai complementaire du banc, s'il y figure
    for s in NSB.get("SCENARIOS_COMPLEMENTAIRES", []):
        if s["code"] == "S10":
            meme = all(np.array_equal(s[k], S10[k]) for k in ("t", "vin", "gx", "dvref")) and s["R0"] == S10["R0"]
            print(f"S10 construit comme scenario_S10.py ; identique a l'essai complementaire du banc : {meme}")

    TACHES = [(c, m, nom, fl, fc) for c in SCENARIOS for m in METHODES for nom, fl, fc in VARIANTES]
    # les plus longues d'abord (S10, adaptatives) pour equilibrer les processus
    TACHES.sort(key=lambda t: (t[0] != "S10", t[1] not in ("PINN-PID", "ELM-PID")))
    print(f"{len(TACHES)} simulations, quatre processus.", flush=True)
    RES = {}
    with mp.get_context("fork").Pool(4) as pool:
        for i, (cle, (iae, d)) in enumerate(pool.imap_unordered(lancer, TACHES, chunksize=1), 1):
            RES[cle] = iae
            print(f"  [{i:3d}/{len(TACHES)}] {cle[0]:4s} {cle[1]:16s} {cle[2]:10s} IAE {iae:10.5f} mV.s "
                  f"({d:5.1f} s)", flush=True)
    print(f"Simulations terminees en {time.time() - T_DEBUT:.0f} s.")
    IAE = {c: {m: {nom: RES[(c, m, nom)] for nom in NOMS_V} for m in METHODES} for c in SCENARIOS}

    # %% Controle prealable : nominal = metriques_banc.csv a 1e-9 relatif
    REF = lire_csv()
    print("\nControle : IAE nominale contre metriques_banc.csv (tolerance 1e-9 relatif)")
    ecart_max, echecs = 0.0, []
    for c in SCENARIOS:
        for m in METHODES:
            r = REF[(c, m)]
            e = abs(IAE[c][m]["nominal"] - r) / r
            ecart_max = max(ecart_max, e)
            if e > 1e-9:
                echecs.append(f"{c} {m} : {IAE[c][m]['nominal']:.10g} contre {r:.10g} (ecart {e:.2e})")
    print(f"  ecart relatif maximal {ecart_max:.2e}")
    if echecs:
        print("  ARRET : le calcul n'est pas valide.\n  " + "\n  ".join(echecs))
        with open(os.path.join(DOSSIER, "separation_resultats.json"), "w", encoding="utf-8") as f:
            json.dump({"valide": False, "echecs_controle": echecs, "IAE_mVs": IAE}, f, indent=1)
        raise SystemExit(1)

    # %% Paires, groupes ou matrice
    SORTIE = {"regle": "signe de la difference d'IAE identique dans les cinq variantes", "valide": True, "ecart_controle_max": ecart_max,
              "variantes": NOMS_V, "IAE_mVs": IAE, "scenarios": {}}
    for c in SCENARIOS:
        ordre = sorted(METHODES, key=lambda m: IAE[c][m]["nominal"])
        print("\n" + "=" * 100 + f"\n {c} : IAE de classement (mV.s) par variante\n" + "=" * 100)
        print("  " + " " * 16 + "".join(f"{n:>12s}" for n in NOMS_V))
        for m in ordre:
            print(f"  {m:16s}" + "".join(f"{IAE[c][m][n]:12.4f}" for n in NOMS_V))
        paires, separe = [], {}
        print(f"\n  Paires (A meilleure au nominal ; ecart = IAE_B / IAE_A - 1, en %)")
        for a, b in itertools.combinations(ordre, 2):
            ec = {n: 100.0 * (IAE[c][b][n] / IAE[c][a][n] - 1.0) for n in NOMS_V}
            sg = [signe(IAE[c][b][n] - IAE[c][a][n]) for n in NOMS_V]
            sep = len(set(sg)) == 1 and sg[0] != 0
            separe[frozenset((a, b))] = sep
            paires.append({"A": a, "B": b, "separee": sep, "ecart_nominal_pct": ec["nominal"],
                           "ecart_min_pct": min(ec.values()), "ecart_max_pct": max(ec.values()), "ecart_pct": ec})
            print(f"    {a:16s} / {b:16s} {'SEPAREE    ' if sep else 'NON SEPAREE'}  nominal {ec['nominal']:+8.3f} %"
                  f"   min {min(ec.values()):+8.3f} %   max {max(ec.values()):+8.3f} %")
        # transitivite de "non separees"
        non_trans = []
        for x, y, z in itertools.permutations(ordre, 3):
            if (not separe[frozenset((x, y))] and not separe[frozenset((y, z))] and separe[frozenset((x, z))]
                    and x < z):
                non_trans.append([x, y, z])
        groupes = None
        if not non_trans:
            groupes, vus = [], set()
            for m in ordre:
                if m in vus:
                    continue
                g = [m] + [n for n in ordre if n != m and not separe[frozenset((m, n))]]
                vus.update(g)
                groupes.append(g)
            # groupes ordonnes seulement s'ils sont contigus dans l'ordre nominal
            plat = [m for g in groupes for m in g]
            if plat != ordre:
                groupes = None
                non_trans.append(["groupes non contigus dans l'ordre nominal"])
        matrice = {a: {b: (None if a == b else separe[frozenset((a, b))]) for b in ordre} for a in ordre}
        SORTIE["scenarios"][c] = {"ordre_nominal": ordre, "paires": paires, "transitive": not non_trans,
                                  "triplets_non_transitifs": non_trans, "groupes_ordonnes": groupes,
                                  "matrice_separee": matrice}
        if groupes is not None:
            print("\n  Groupes ordonnes : " + " < ".join("{" + ", ".join(g) + "}" for g in groupes))
        else:
            print("\n  Relation non transitive (" + "; ".join(" ~ ".join(t) for t in non_trans) + ")")
            print("  Matrice des paires (S = separee, - = non separee) :")
            print("  " + " " * 16 + "".join(f"{m[:10]:>11s}" for m in ordre))
            for a in ordre:
                print(f"  {a:16s}" + "".join(f"{'':>11s}" if a == b else f"{('S' if matrice[a][b] else '-'):>11s}"
                                            for b in ordre))

    with open(os.path.join(DOSSIER, "separation_resultats.json"), "w", encoding="utf-8") as f:
        json.dump(SORTIE, f, indent=1)

    print("\n" + "=" * 100 + "\n Resume\n" + "=" * 100)
    for c in SCENARIOS:
        r = SORTIE["scenarios"][c]
        if r["groupes_ordonnes"] is not None:
            print(f"  {c:4s} : " + " < ".join("{" + ", ".join(g) + "}" for g in r["groupes_ordonnes"]))
        else:
            ns_ = [f"{p['A']} ~ {p['B']}" for p in r["paires"] if not p["separee"]]
            print(f"  {c:4s} : non transitive ; paires non separees : " + "; ".join(ns_))
    print(f"\nseparation_resultats.json ecrit. Duree totale {time.time() - T_DEBUT:.0f} s.")
