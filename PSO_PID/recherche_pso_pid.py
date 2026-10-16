# =============================================================================
# recherche_pso_pid.py
#
# OBJECTIF
#   Regler hors ligne les gains P, I, D du bloc "PID Controller" de
#   Ziegler-Nichols par essaim de particules (PSO). L'essaim est celui de
#   Gaing (2004) ; le cout, le traitement des contraintes et les conditions
#   initiales du reglage sont les suivants, chacun appuye sur la
#   litterature :
#     cout  IAE (critere cite par Gaing, eq. 6) sur quatre essais de reglage
#           qui contiennent des perturbations de charge, de Vin, de consigne
#           et du bruit de mesure (Astrom, Panagopoulos, Hagglund 1998 ;
#           Panagopoulos, Astrom, Hagglund 2002 ; Krohling, Rey 2001 ;
#           Alipoor 2009) ; chaque IAE est rapportee a celle de
#           Ziegler-Nichols sur le meme essai ; le demarrage pese 3/13 et le
#           reste 10/13 ;
#     contraintes  marge de phase >= 30 degres, coupure <= fs/10 aux 12
#           coins, traitees par les regles de faisabilite de Deb (2000),
#           appliquees au PSO comme Toscano Pulido et Coello Coello (2004) ;
#     conditions initiales  pendant le reglage, le PID part des conditions
#           initiales du bloc (integrateur 0.5, filtre 0.01), celles du
#           fonctionnement reel.
#   Essaim (Gaing 2004) : population 50, 100 iterations, inertie de 0.9 a
#   0.4 (eq. 12), c1 = c2 = 2, |v| <= 2 (moitie de la plage), positions
#   ramenees dans [0 ; 4] fois les gains de Ziegler-Nichols (eq. 15).
#   Cinq graines. Les onze essais de developpement ne servent qu'a juger
#   (banc_pso_pid.py).
#
# DEUX FACONS DE LE LANCER
#   RECALCULER = False (reglage fourni) : reprend recherche_pso_pid.json et
#   recalcule le cout de chaque gbest pour controle (une minute).
#   RECALCULER = True : refait les cinq recherches (35 a 70 minutes selon la
#   machine). Sur une autre machine, les arrondis peuvent faire diverger
#   l'essaim apres quelques iterations : le gbest peut differer un peu.
#
# CE QUE PRODUIT CE SCRIPT
#   recherche_pso_pid.json : pour chaque graine, gbest, cout, rapports a
#   Ziegler-Nichols sur les essais de reglage, historique de convergence.
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT
#   python recherche_pso_pid.py
#
# ORDRE D'EXECUTION
#   1. ce script ; 2. banc_pso_pid.py ; 3. Construction_PSO_PID.m ;
#   4. Simuler_PSO_PID.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques, reglages et definitions du banc commun

import os                                     # chemins de fichiers
import json                                   # lecture et ecriture au format texte
import time                                   # duree d'execution
import multiprocessing as mp                  # evaluation en parallele des individus
import numpy as np                            # calcul numerique
from scipy.signal import cont2discrete        # modele moyen discretise (marges)

RECALCULER = False                            # True : refaire les cinq recherches

try:                                          # dossier du script
    DOSSIER_PSO = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PSO = os.getcwd()
with open(os.path.join(DOSSIER_PSO, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC or "iL_max_dem_A" not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py ne contient pas les definitions attendues (marque de fin, observation du courant).")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])  # constantes, Circuit, simuler, PIDClassique

FICHIER = os.path.join(DOSSIER_PSO, "recherche_pso_pid.json")

# Essaim (Gaing 2004)
N_POP, N_ITER = 50, 100                       # population, iterations
W_MAX, W_MIN = 0.9, 0.4                       # inertie (eq. 12)
C1 = C2 = 2.0                                 # constantes d'acceleration
K_MIN, K_MAX = 0.0, 4.0                       # bornes des multiplicateurs de ZN 
V_MAX = K_MAX / 2.0                           # vitesse maximale (moitie de la plage)
GRAINES = (1, 2, 3, 4, 5)
# Contraintes (base commune : celles du meilleur PID fige)
MARGE_MINI = 30.0                             # degres
FC_MAXI = FSW / 10.0                          # Hz
# Cout
POIDS_DEM = 3.0 / 13.0                        # poids du demarrage (3 termes sur 13 dans J)
K30 = int(round(0.030 / TC))                  # fin du demarrage (comme J)
GRAINE_BRUIT_REGLAGE = 7                      # bruit du reglage : autre graine que S7b (20261002)
FC_BRUIT, ECART_BRUIT = 10e3, 0.01            # meme nature que S7b
PAS_CAN = 2.5 / 4096 * 80                     # meme CAN que S7a (48.8 mV)


class PIDParallele:
    """Bloc "PID Controller" de Simulink (forme Parallel, Forward Euler,
    clamping, sortie [0.01 ; 0.99]) avec des gains P, I, D donnes ; memes
    operations que PIDClassique du banc commun ; conditions initiales du
    bloc."""

    def __init__(self, P, I, D):
        self.P, self.I, self.D = P, I, D
        self.xI, self.xF = PID_CI_INTEGRATEUR, PID_CI_FILTRE

    def pas(self, e):
        derivee = PID_N * (self.D * e - self.xF)
        u = self.P * e + self.xI + derivee
        u_sat = min(max(u, D_MIN), D_MAX)
        entree_I = self.I * e
        if not ((u != u_sat) and (np.sign(entree_I) == np.sign(u - u_sat))):
            self.xI += TC * entree_I
        self.xF += TC * derivee
        return u_sat


# %% ETAPE 1 : contraintes (marge et coupure aux 12 coins, comme recherche_pid_fige.py)

W_PULS = np.logspace(np.log10(2 * np.pi * 20.0), np.log10(np.pi / TC * 0.999), 60000)   # pulsations (rad/s)
Z_CERCLE = np.exp(1j * W_PULS * TC)
COINS = [(R, V) for R in (4.0, 5.0, 25.0, 98.0) for V in (160.0, 200.0, 240.0)]
REPONSES = {}
for R, V in COINS:
    num_d, den_d, _ = cont2discrete(([V], [L_BOB * C_CONV, L_BOB / R, 1.0]), TC, method="zoh")
    REPONSES[(R, V)] = np.polyval(np.squeeze(num_d), Z_CERCLE) / np.polyval(den_d, Z_CERCLE)


def frequentiel(a, b, c):
    """Plus petite marge de phase (degres) et plus haute frequence de coupure
    (Hz) aux 12 coins (meme calcul que recherche_pid_fige.py)."""
    P, I, D = PID_P * a, PID_I * b, PID_D * c
    Cz = P + I * TC / (Z_CERCLE - 1) + D * PID_N / (1 + PID_N * TC / (Z_CERCLE - 1))
    pire, fc = 180.0, 0.0
    for G in REPONSES.values():
        boucle = Cz * G
        module = np.abs(boucle)
        phase = np.degrees(np.unwrap(np.angle(boucle)))
        croisements = np.where(np.diff(np.sign(module - 1.0)))[0]
        for j in croisements:
            pire = min(pire, (phase[j] + 360.0) % 360.0 - 180.0)
            fc = max(fc, W_PULS[j] / (2 * np.pi))
    return float(pire), float(fc)



# %% ETAPE 2 : les quatre essais de reglage et le cout (aucun n'est un des onze essais)

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
    if bruit:
        blanc = np.random.default_rng(GRAINE_BRUIT_REGLAGE).standard_normal(n)
        a = np.exp(-2.0 * np.pi * FC_BRUIT * TC)
        y, etat = np.zeros(n), 0.0
        for k in range(n):
            etat = a * etat + (1.0 - a) * blanc[k]
            y[k] = etat
        sc["bruit"] = ECART_BRUIT * y / np.std(y)
    return sc


REGLAGE = [
    # E1 : 7 ohms, 190 V ; charge portee a 4.5 ohms sur [40 ; 65) ms
    essai_reglage("E1", 7.0, 190.0, 0.090, (0.040, 0.065), gx=[(0.040, 0.065, 1 / 4.5 - 1 / 7.0)]),
    # E2 : 15 ohms, 210 V ; Vin 175 V sur [40 ; 65) ms puis 235 V
    essai_reglage("E2", 15.0, 210.0, 0.090, (0.040, 0.065), vin=[(0.040, 0.065, 175.0), (0.065, 1.0, 235.0)]),
    # E3 : 60 ohms, 220 V ; echelon de consigne +10 V sur [40 ; 65) ms
    essai_reglage("E3", 60.0, 220.0, 0.090, (0.040, 0.065), dvref=[(0.040, 0.065, 10.0)]),
    # E4 : 9 ohms, 200 V ; bruit (10 mV, 10 kHz, autre graine) et CAN 12 bits ensemble
    essai_reglage("E4", 9.0, 200.0, 0.060, (), bruit=True, q=PAS_CAN),
]
for sc in REGLAGE:                         # aucun point de reglage n'est un point des onze essais
    for s in SCENARIOS:
        if abs(s["R0"] - sc["R0"]) < 1e-9 and np.any(np.abs(s["vin"] - sc["vin"][0]) < 1e-9):
            raise RuntimeError(f"Le point de reglage {sc['code']} apparait dans l'essai {s['code']}.")


def iae(P, I, D, sc):
    """IAE du demarrage (0 a 30 ms) et IAE de 30 ms a la fin, conditions
    initiales du bloc."""
    sim = simuler(PIDParallele(P, I, D), sc)
    e = np.abs(sim["consigne"] - sim["v"])
    return float(np.sum(e[:K30]) * TC), float(np.sum(e[K30:]) * TC)


IAE_ZN = None                                 # normalisation (calculee une fois, dans chaque processus)


def iae_zn():
    global IAE_ZN
    if IAE_ZN is None:
        IAE_ZN = np.array([iae(PID_P, PID_I, PID_D, sc) for sc in REGLAGE])
    return IAE_ZN


def violation(marge, fc):
    """Violation des contraintes (0 si admissible), normalisee (Deb 2000)."""
    return max(0.0, (MARGE_MINI - marge) / MARGE_MINI) + max(0.0, (fc - FC_MAXI) / FC_MAXI)


def evaluer(x):
    """Contraintes, puis, si admissible, le cout J_reglage."""
    a, b, c = (float(t) for t in x)
    marge, fc = frequentiel(a, b, c)
    viol = violation(marge, fc)
    res = {"a": a, "b": b, "c": c, "marge_min": marge, "fc_max": fc, "violation": viol, "admissible": viol == 0.0}
    if res["admissible"]:
        r = np.array([iae(PID_P * a, PID_I * b, PID_D * c, sc) for sc in REGLAGE]) / iae_zn()
        res["ratios"] = r.tolist()
        res["J"] = float(POIDS_DEM * np.mean(r[:, 0]) + (1.0 - POIDS_DEM) * np.mean(r[:, 1]))
    return res


def meilleur(r1, r2):
    """Regles de Deb : r1 bat-il r2 ? Admissible bat non admissible ; deux
    admissibles : plus petit J ; deux non admissibles : plus petite violation."""
    if r1["admissible"] and r2["admissible"]:
        return r1["J"] < r2["J"]
    if r1["admissible"] != r2["admissible"]:
        return r1["admissible"]
    return r1["violation"] < r2["violation"]


# %% ETAPE 3 : l'essaim (Gaing, etapes 1 a 9, comparaisons par les regles de Deb)

def recherche(graine, pool):
    rng = np.random.default_rng(graine)
    x = rng.uniform(K_MIN, K_MAX, (N_POP, 3))
    v = rng.uniform(-V_MAX, V_MAX, (N_POP, 3))
    res = pool.map(evaluer, [tuple(r) for r in x])
    pbest, res_pbest = x.copy(), list(res)
    g = 0
    for i in range(1, N_POP):
        if meilleur(res_pbest[i], res_pbest[g]):
            g = i

    def note(it):
        adm = [r["J"] for r in res if r["admissible"]]
        return {"iter": it, "J_gbest": res_pbest[g].get("J"), "violation_gbest": res_pbest[g]["violation"],
                "gbest": pbest[g].tolist(), "part_admissible": len(adm) / N_POP,
                "J_med_admissibles": float(np.median(adm)) if adm else None}

    histo = [note(0)]
    for it in range(1, N_ITER + 1):
        w = W_MAX - (W_MAX - W_MIN) / N_ITER * it
        r1 = rng.random((N_POP, 3))
        r2 = rng.random((N_POP, 3))
        v = w * v + C1 * r1 * (pbest - x) + C2 * r2 * (pbest[g] - x)
        v = np.clip(v, -V_MAX, V_MAX)
        x = np.clip(x + v, K_MIN, K_MAX)
        res = pool.map(evaluer, [tuple(r) for r in x])
        for i in range(N_POP):
            if meilleur(res[i], res_pbest[i]):
                pbest[i], res_pbest[i] = x[i].copy(), res[i]
        for i in range(N_POP):
            if meilleur(res_pbest[i], res_pbest[g]):
                g = i
        histo.append(note(it))
    rg = res_pbest[g]
    return {"graine": graine, "gbest": pbest[g].tolist(), "J_reglage": rg.get("J"), "admissible": rg["admissible"],
            "gains": {"P": PID_P * pbest[g][0], "I": PID_I * pbest[g][1], "D": PID_D * pbest[g][2]},
            "marge_min": rg["marge_min"], "fc_max": rg["fc_max"],
            "ratios": {sc["code"]: {"demarrage": q[0], "apres_30ms": q[1]}
                       for sc, q in zip(REGLAGE, rg.get("ratios", []))},
            "historique": histo}


# %% ETAPE 4 : les cinq recherches (ou leur controle)

if __name__ == "__main__":
    T0 = time.time()
    print("Essais de reglage : " + " ; ".join(f"{sc['code']} R = {sc['R0']:g} ohms, Vin = {sc['vin'][0]:g} V, "
                                               f"{sc['duree'] * 1e3:.0f} ms" for sc in REGLAGE))
    print("IAE de Ziegler-Nichols (mV.s, demarrage / apres 30 ms) : " +
          " ; ".join(f"{sc['code']} {q[0] * 1e3:.2f} / {q[1] * 1e3:.2f}" for sc, q in zip(REGLAGE, iae_zn())))
    ZN = evaluer((1.0, 1.0, 1.0))
    print(f"Ziegler-Nichols : marge minimale {ZN['marge_min']:.1f} degres, coupure {ZN['fc_max']:.0f} Hz, "
          f"{'admissible' if ZN['admissible'] else 'refuse par les contraintes'} (il ne sert qu'a normaliser).")
    resultats = {}
    if os.path.isfile(FICHIER) and not RECALCULER:
        with open(FICHIER, encoding="utf-8") as f:
            resultats = json.load(f)["recherches"]
        print(f"{len(resultats)} recherche(s) reprise(s) de recherche_pso_pid.json.")
    with mp.Pool(min(4, os.cpu_count() or 1)) as pool:
        for graine in GRAINES:
            cle = f"graine={graine}"
            if cle in resultats:
                r = resultats[cle]
                J_recalc = evaluer(tuple(r["gbest"]))["J"]
                ecart = abs(J_recalc - r["J_reglage"]) / r["J_reglage"]
                print(f"  {cle} : repris, J_reglage {r['J_reglage']:.6f}, recalcule {J_recalc:.6f} (ecart {ecart:.1e})")
                if ecart > 1e-9:
                    raise RuntimeError(f"{cle} : le cout recalcule differe du fichier.")
                continue
            t = time.time()
            r = recherche(graine, pool)
            resultats[cle] = r
            with open(FICHIER, "w", encoding="utf-8") as f:
                json.dump({"format": 1, "methode": "PSO-PID : essaim de Gaing 2004",
                           "recherches": resultats}, f, indent=1, allow_nan=False)
            print(f"  {cle} : J_reglage = {r['J_reglage']:.5f}, gains x ZN ({r['gbest'][0]:.4f}, {r['gbest'][1]:.4f}, "
                  f"{r['gbest'][2]:.4f}), marge {r['marge_min']:.1f} degres, coupure {r['fc_max']:.0f} Hz "
                  f"({time.time() - t:.0f} s)", flush=True)
    print(f"Duree totale : {time.time() - T0:.0f} s.")
