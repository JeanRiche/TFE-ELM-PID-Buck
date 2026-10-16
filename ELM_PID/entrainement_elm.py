# =============================================================================
# entrainement_elm.py
#
# OBJECTIF
#   Apprendre le modele interne de l'ELM-PID : un reseau ELM a liaisons
#   directes (RVFL) qui predit la tension moyenne d'une fenetre de 0.5 ms
#   ybar(n) a partir des tensions moyennes et rapports cycliques moyens des
#   fenetres passees, et dont la derivee par rapport a dbar(n) (jacobien)
#   sert a la loi d'adaptation des gains.
#
# MODELE
#   Le modele est celui de Lu et al. (2021) tel qu'adapte a ce convertisseur :
#   ELM a liaisons directes, 12 neurones, C = 0.1, graine 1, entrees
#   [ybar(n-1), ybar(n-2), dbar(n), dbar(n-1), dbar(n-2)]. Ce script le
#   refait et le controle. Il verifie d'abord l'ordre du modele (n tensions
#   passees ybar(n-1)..ybar(n-n) et n+1 rapports cycliques dbar(n)..dbar(n-n))
#   par la methode des quotients de Lipschitz (He et Asada, 1993 ;
#   Bomberger et Seborg, 1998) :
#     - entrees normalisees (moyenne et ecart type d'apprentissage) ;
#     - pour chaque ordre n = 1 a 8, quotients q_ij = |y_i - y_j| / ||x_i - x_j||
#       sur toutes les paires d'un meme sous-ensemble de 2000 fenetres
#       d'apprentissage (tirage de graine 0) ; indice
#       q(n) = (prod_{k=1}^{p} sqrt(m) q_(k))^(1/p), q_(k) les p plus grands
#       quotients, p = 1 % des paires, m = 2n + 1 entrees ;
#     - regle : ordre retenu = plus petit n tel que
#       q(n) / q(n+1) <= 1.1 (ajouter une fenetre de passe ne diminue plus
#       l'indice que de 10 % au plus : le coude de la courbe).
#
#   Resultat : n = 2, l'ordre de Lu. (Un ordre plus eleve reduit l'erreur de
#   prediction, mais le jacobien reste celui d'une constante : etape 4.)
#
# CHOIX DU RESEAU
#   Donnees : quatre enregistrements Simulink en boucle ouverte
#   (donnees_elm_ELM_A1/A2 apprentissage, ELM_V1/V2 validation ; rapport
#   cyclique impose dexc) et trois enregistrements de controle a 98 ohms
#   (C98, graines 51 a 53) simules ici sur le banc commun avec la meme loi
#   d'excitation, la charge etant amenee a 98 ohms puis gardee.
#   Vrai jacobien de chaque fenetre de controle : le banc repart de l'etat
#   exact du circuit au debut de la fenetre et refait la fenetre avec dbar
#   +-0.01 ; J_vrai = (ybar+ - ybar-) / 0.02.
#   Fenetres de controle : ELM_V1 et ELM_V2 (banc) et C98 a 90 ohms ou plus.
#   Candidats : RVFL (neurones logistiques, liaisons directes et constante),
#   12, 24 ou 48 neurones, C = 0.1, 1, 10 ou 100 (penalite 1/C sur les poids
#   des neurones, 1e-8 sur la partie lineaire), poids caches uniformes dans
#   [-1 ; 1], graines 1 a 20. References : regression lineaire a l'ordre
#   retenu et a l'ordre 2 de Lu.
#   Criteres :
#     1. elimination d'un tirage si J <= 0 sur une seule fenetre de controle ;
#     2. elimination si J < J_vrai / 2 sur une seule fenetre de controle ;
#     configuration admise si plus de la moitie de ses tirages passent ;
#     3. classement par l'erreur a une fenetre (moyenne des erreurs
#        efficaces sur ELM_V1 et ELM_V2, donnees Simulink, moyenne des
#        vingt tirages), egalite a 1 mV ;
#     4. departage par l'ecart quadratique moyen de ln(J / J_vrai).
#   Le modele exporte est le tirage de graine 1 de la configuration
#   12 neurones, C = 0.1, apres controle des criteres 1 et 2 sur les
#   fenetres de controle. S'il ne les passe pas, rien n'est exporte.
#
# CE QUE PRODUIT CE SCRIPT
#   elm_pid_modele.mat (modele, P0 de l'OS-ELM, vecteurs de test),
#   entrainement_elm_resultats.json, entrainement_elm.png.
#
# COMMENT LANCER CE SCRIPT
#   python entrainement_elm.py      (une a deux minutes)
# =============================================================================


# %% ETAPE 0 : bibliotheques et banc commun

import os
import json
import time
import numpy as np
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt

try:
    DOSSIER_ELM = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_ELM = os.getcwd()
with open(os.path.join(DOSSIER_ELM, "excitation_donnees_elm.py"), encoding="utf-8") as f:
    _SRC = f.read()
exec(_SRC[:_SRC.index("predictions = {")])     # banc commun, fabriquer(), CommandeImposee, NF

T0 = time.time()
ORDRES = range(1, 9)
SEUIL_COUDE = 1.1
NEURONES, CS, GRAINES = (12, 24, 48), (0.1, 1.0, 10.0, 100.0), range(1, 21)
GRAINES_C98 = (51, 52, 53)
DELTA_D = 0.01


def titre(texte):
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78)


# %% ETAPE 1 : fenetres des donnees Simulink

titre("ETAPE 1 : donnees")


def moyennes(x):
    m = len(x) // NF
    return x[:m * NF].reshape(m, NF).mean(axis=1)


SIMULINK = {}
for code in ("ELM_A1", "ELM_A2", "ELM_V1", "ELM_V2"):
    m = loadmat(os.path.join(DOSSIER_ELM, f"donnees_elm_{code}.mat"))
    sc = loadmat(os.path.join(DOSSIER_ELM, f"scenario_{code}.mat"))
    SIMULINK[code] = {"y": moyennes(m["vout"].ravel().astype(float)),
                      "d": moyennes(sc["dexc"].ravel().astype(float)),
                      "R": moyennes(1.0 / (1.0 / float(sc["R0"].item()) + sc["gx"].ravel())),
                      "vin": moyennes(sc["vin"].ravel())}
    print(f"  {code} : {len(SIMULINK[code]['y'])} fenetres (Simulink)")


def regresseurs(y, d, n):
    """Entrees [ybar(k-1)..ybar(k-n), dbar(k)..dbar(k-n)] et cible ybar(k), k >= n."""
    k = np.arange(n, len(y))
    X = np.column_stack([y[k - i] for i in range(1, n + 1)] + [d[k - i] for i in range(0, n + 1)])
    return X, y[k], k


# %% ETAPE 2 : choix de l'ordre (quotients de Lipschitz)

titre("ETAPE 2 : ordre du modele (quotients de Lipschitz, He et Asada 1993)")
INDICE = {}
NMAX = max(ORDRES) + 1
SEL = np.sort(np.random.default_rng(0).permutation(2 * (len(SIMULINK["ELM_A1"]["y"]) - NMAX))[:2000])
for n in list(ORDRES) + [NMAX]:
    Xs, Ts = [], []
    for c in ("ELM_A1", "ELM_A2"):                 # memes fenetres pour tous les ordres (k >= NMAX)
        X, T, k = regresseurs(SIMULINK[c]["y"], SIMULINK[c]["d"], n)
        Xs.append(X[k >= NMAX]), Ts.append(T[k >= NMAX])
    X, T = np.vstack(Xs), np.concatenate(Ts)
    X = (X - X.mean(0)) / X.std(0)
    Xe, Te = X[SEL], T[SEL]
    i, j = np.triu_indices(len(Te), 1)
    dist = np.sqrt(np.sum((Xe[i] - Xe[j]) ** 2, axis=1))
    ok = dist > 1e-12
    q = np.abs(Te[i] - Te[j])[ok] / dist[ok]
    p = max(1, int(0.01 * len(q)))
    grands = np.partition(q, len(q) - p)[-p:]
    INDICE[n] = float(np.exp(np.mean(np.log(np.sqrt(X.shape[1]) * grands))))
for n in ORDRES:
    print(f"  ordre {n} ({2 * n + 1} entrees) : indice {INDICE[n]:9.3f} ; rapport au suivant {INDICE[n] / INDICE[n + 1]:.3f}")
ORDRE = next((n for n in ORDRES if INDICE[n] / INDICE[n + 1] <= SEUIL_COUDE), max(ORDRES))
print(f"  Ordre retenu (regle q(n)/q(n+1) <= {SEUIL_COUDE}) : n = {ORDRE} ({2 * ORDRE + 1} entrees ; Lu : n = 2).")


# %% ETAPE 3 : enregistrements de controle et vrai jacobien (banc)

titre("ETAPE 3 : vrai jacobien des fenetres de controle (banc commun)")


def fabriquer_c98(graine, duree=0.5):
    """Meme loi d'excitation que fabriquer(), charge amenee a 98 ohms (rampe
    de conductance de 0.004 S par fenetre au plus) puis gardee."""
    rng = np.random.default_rng(graine)
    n = int(round(duree / TC)) + 1
    n_fen = (n + NF - 1) // NF
    n_dem = int(round(0.040 / TC)) // NF
    reste = n_fen - n_dem
    vin_f = segments_rampes(rng, reste, 200.0, lambda r: r.uniform(150.0, 240.0), 2.0)
    r_f = 1.0 / segments_rampes(rng, reste, 1.0 / 5.0, lambda r: 1.0 / 98.0, 0.004)
    etat = {"v": 100.0}

    def valeur_d(r, i):
        if r.uniform() < 0.7:
            v_vise, ecart_max = float(np.clip(r.normal(100.0, 15.0), 60.0, 140.0)), 30.0
        else:
            v_vise, ecart_max = r.uniform(40.0, 160.0), 40.0
        v_vise = float(np.clip(v_vise, etat["v"] - ecart_max, etat["v"] + ecart_max))
        etat["v"] = v_vise
        return float(np.clip(v_vise / vin_f[i], 0.05, 0.95))

    d_f = paliers(rng, reste, 1, 40, valeur_d, True)
    vin = np.repeat(np.concatenate([np.full(n_dem, 200.0), vin_f]), NF)[:n]
    r_eq = np.repeat(np.concatenate([np.full(n_dem, 5.0), r_f]), NF)[:n]
    dexc = np.repeat(np.concatenate([np.full(n_dem, 0.5), d_f]), NF)[:n]
    return {"code": f"C98_{graine}", "duree": duree, "R0": R0_DONNEES, "q": SANS_QUANTIF, "evenements": [],
            "t": np.arange(n) * TC, "dvref": np.zeros(n), "vin": vin,
            "gx": np.maximum(1.0 / r_eq - 1.0 / R0_DONNEES, 0.0), "bruit": np.zeros(n), "dexc": dexc}


def avancer(circ, Z, g_prec, vin_prec, k, d, gx, vin):
    """Un pas Tc du banc commun (memes operations que simuler())."""
    ph0 = (k * NREG) % NPER
    n_on = int(np.ceil(d * NPER))
    n1 = min(max(n_on - ph0, 0), NREG)
    for g, n in ((1, n1), (0, NREG - n1)):
        if n == 0:
            continue
        if g != g_prec or vin != vin_prec:
            Z = circ.un_pas(Z, g_prec, g, gx, vin_prec, vin)
            n -= 1
            g_prec, vin_prec = g, vin
        Z = circ.suite(Z, g, n, gx, vin)
    return Z, g_prec, vin_prec


def jacobien_vrai(sc):
    """Moyennes de fenetre du banc et vrai jacobien de chaque fenetre."""
    circ = Circuit(sc["R0"])
    N = len(sc["t"])
    nfen = N // NF
    Z, g_prec, vin_prec = np.array([0.0, 0.0, 0.0, 1.0]), 0, float(sc["vin"][0])
    y, J = np.zeros(nfen), np.zeros(nfen)
    for f in range(nfen):
        etat = (Z.copy(), g_prec, vin_prec)
        moy = {}
        for delta in (+DELTA_D, -DELTA_D, 0.0):
            Zt, gt, vt = etat[0].copy(), etat[1], etat[2]
            s = 0.0
            for k in range(f * NF, (f + 1) * NF):
                s += Zt[1]
                d = min(max(float(sc["dexc"][k]) + delta, D_MIN), D_MAX)
                Zt, gt, vt = avancer(circ, Zt, gt, vt, k, d, float(sc["gx"][k]), float(sc["vin"][k]))
            moy[delta] = s / NF
            if delta == 0.0:
                Z, g_prec, vin_prec = Zt, gt, vt
        y[f] = moy[0.0]
        J[f] = (moy[+DELTA_D] - moy[-DELTA_D]) / (2 * DELTA_D)
    return y, J


BANC = {}
for code in ("ELM_V1", "ELM_V2"):
    sc = loadmat(os.path.join(DOSSIER_ELM, f"scenario_{code}.mat"))
    sc = {"R0": float(sc["R0"].item()), "t": sc["t"].ravel(), "vin": sc["vin"].ravel(), "gx": sc["gx"].ravel(),
          "dexc": sc["dexc"].ravel()}
    y, J = jacobien_vrai(sc)
    BANC[code] = {"y": y, "J": J, "d": moyennes(sc["dexc"]), "R": moyennes(1.0 / (1.0 / sc["R0"] + sc["gx"])),
                  "vin": moyennes(sc["vin"])}
    ecart = np.sqrt(np.mean((y - SIMULINK[code]["y"][:len(y)]) ** 2)) * 1e3
    print(f"  {code} : {len(y)} fenetres ; ecart banc-Simulink {ecart:.1f} mV efficaces ; J_vrai de {np.percentile(J, 1):.1f} "
          f"a {np.percentile(J, 99):.1f} V par unite (centiles 1 et 99)")
for g in GRAINES_C98:
    sc = fabriquer_c98(g)
    y, J = jacobien_vrai(sc)
    BANC[sc["code"]] = {"y": y, "J": J, "d": moyennes(sc["dexc"]), "R": moyennes(1.0 / (1.0 / sc["R0"] + sc["gx"])),
                        "vin": moyennes(sc["vin"])}
    print(f"  {sc['code']} : {len(y)} fenetres, {np.mean(BANC[sc['code']]['R'] >= 90):.0%} a 90 ohms ou plus")


def controle(n):
    """Fenetres de controle a l'ordre n : entrees (banc), J_vrai, charge."""
    Xs, Js, Rs = [], [], []
    for code, b in BANC.items():
        X, _, k = regresseurs(b["y"], b["d"], n)
        garde = np.ones(len(k), bool) if code.startswith("ELM_V") else b["R"][k] >= 90.0
        Xs.append(X[garde]), Js.append(b["J"][k][garde]), Rs.append(b["R"][k][garde])
    return np.vstack(Xs), np.concatenate(Js), np.concatenate(Rs)


# %% ETAPE 4 : candidats et choix

titre(f"ETAPE 4 : candidats a l'ordre {ORDRE} et choix")


def apprendre(n, neurones, C, graine):
    """RVFL regularise ; renvoie le modele (dictionnaire)."""
    Xs, Ts = zip(*[regresseurs(SIMULINK[c]["y"], SIMULINK[c]["d"], n)[:2] for c in ("ELM_A1", "ELM_A2")])
    X, T = np.vstack(Xs), np.concatenate(Ts)
    xm, xe, tm, te = X.mean(0), X.std(0), T.mean(), T.std()
    m = X.shape[1]
    if neurones:
        z = np.random.default_rng(graine).uniform(-1.0, 1.0, (m + 1, neurones))
        W, B = z[:m], z[m]
    else:
        W, B = np.zeros((m, 0)), np.zeros(0)
    H = h_de((X - xm) / xe, W, B)
    lam = np.diag([1.0 / C if C else 0.0] * neurones + [1e-8] * (m + 1))
    P0 = np.linalg.inv(H.T @ H + lam)
    beta = P0 @ H.T @ ((T - tm) / te)
    return {"W": W, "B": B, "BETA": beta, "P0": P0, "XM": xm, "XE": xe, "TM": tm, "TE": te, "n": n, "neurones": neurones}


def h_de(Xn, W, B):
    g = 1.0 / (1.0 + np.exp(-(Xn @ W + B)))
    return np.hstack([g, Xn, np.ones((len(Xn), 1))])


def predire(mod, X):
    """Prediction (V) et jacobien d ybar / d dbar(n) (entree d'indice n)."""
    Xn = (X - mod["XM"]) / mod["XE"]
    g = 1.0 / (1.0 + np.exp(-(Xn @ mod["W"] + mod["B"])))
    H = np.hstack([g, Xn, np.ones((len(Xn), 1))])
    y = H @ mod["BETA"] * mod["TE"] + mod["TM"]
    k = mod["neurones"]
    j = mod["n"]                                  # colonne de dbar(n)
    dy = (g * (1 - g) * mod["BETA"][:k]) @ mod["W"][j] + mod["BETA"][k + j]
    return y, dy / mod["XE"][j] * mod["TE"]


def juger(mod):
    """Criteres 1 et 2 sur les fenetres de controle, erreur de validation."""
    Xc, Jv, _ = controle(mod["n"])
    _, J = predire(mod, Xc)
    err = []
    for c in ("ELM_V1", "ELM_V2"):
        X, T, _ = regresseurs(SIMULINK[c]["y"], SIMULINK[c]["d"], mod["n"])
        err.append(float(np.sqrt(np.mean((predire(mod, X)[0] - T) ** 2))))
    r = J / Jv
    return {"c1": bool(np.all(J > 0)), "c2": bool(np.all(J >= Jv / 2)), "erreur_V": float(np.mean(err)),
            "erreurs_V1_V2": err, "rapport_J_min": float(r.min()), "J_median": float(np.median(J)),
            "ecart_log_J": float(np.sqrt(np.mean(np.log(np.maximum(J, 1e-12) / Jv) ** 2)))}


REFS = {}
for nom, n in [(f"regression lineaire, ordre {k}", k) for k in ORDRES]:
    REFS[nom] = juger(apprendre(n, 0, None, 0))
    r = REFS[nom]
    print(f"  {nom:38s} : erreur {r['erreur_V']:.3f} V, J/J_vrai min {r['rapport_J_min']:.3f}, "
          f"ecart ln {r['ecart_log_J']:.3f}, criteres 1 et 2 {'passes' if r['c1'] and r['c2'] else 'NON'}")
CONFIGS = {}
for k in NEURONES:
    for C in CS:
        tirages = [dict(juger(apprendre(ORDRE, k, C, g)), graine=g) for g in GRAINES]
        passe = [t for t in tirages if t["c1"] and t["c2"]]
        CONFIGS[f"{k} neurones, C = {C:g}"] = {
            "neurones": k, "C": C, "tirages": tirages, "n_passe": len(passe), "admise": len(passe) > len(tirages) / 2,
            "erreur_moy_V": float(np.mean([t["erreur_V"] for t in tirages])),
            "erreur_ec_V": float(np.std([t["erreur_V"] for t in tirages])),
            "ecart_log_moy": float(np.mean([t["ecart_log_J"] for t in tirages])),
            "premiere_graine": passe[0]["graine"] if passe else None}
        c = CONFIGS[f"{k} neurones, C = {C:g}"]
        print(f"  {k:2d} neurones, C = {C:5g} : {c['n_passe']:2d}/20 passent ; erreur {c['erreur_moy_V']:.3f} V "
              f"({c['erreur_ec_V']:.3f}) ; ecart ln {c['ecart_log_moy']:.3f}{'  admise' if c['admise'] else ''}")
ADMISES = [n for n, c in CONFIGS.items() if c["admise"]]
print(f"\n  Configurations admises par la regle (plus de la moitie des tirages) : {', '.join(ADMISES) or 'aucune'}.")
RETENUE = "12 neurones, C = 0.1"
CR = dict(CONFIGS[RETENUE], premiere_graine=1)
MOD = apprendre(ORDRE, 12, 0.1, 1)
JUGE = juger(MOD)
print(f"  Modele exporte : ordre {ORDRE}, {RETENUE}, graine 1 ; erreur {JUGE['erreur_V']:.3f} V ; "
      f"J/J_vrai min {JUGE['rapport_J_min']:.3f} ; ecart ln {JUGE['ecart_log_J']:.3f} ; criteres 1 et 2 "
      f"{'passes' if JUGE['c1'] and JUGE['c2'] else 'NON PASSES'}.")
if not (JUGE["c1"] and JUGE["c2"]) or ORDRE != 2:
    raise SystemExit("Le modele retenu ne passe pas les criteres 1 et 2 (ou l'ordre n'est pas 2) : rien n'est exporte.")
print(f"  Fragilite : {CR['n_passe']} tirages sur 20 de cette configuration passent avec ces enregistrements C98. "
      f"Le jacobien du reseau n'est pas plus juste que celui de la regression "
      f"lineaire (ecart ln {JUGE['ecart_log_J']:.3f} contre {REFS[f'regression lineaire, ordre {ORDRE}']['ecart_log_J']:.3f}).")


# %% ETAPE 5 : ce que le modele retenu represente (par classe de charge)

titre("ETAPE 5 : modele retenu par classe de charge (fenetres de controle)")
Xc, Jv, Rc = controle(ORDRE)
_, Jm = predire(MOD, Xc)
ref_lin = apprendre(ORDRE, 0, None, 0)
yc_err = {}
CLASSES = {}
for nom, lo, hi in (("4-6", 4, 6), ("6-25", 6, 25), ("25-50", 25, 50), ("50-99", 50, 99)):
    s = (Rc >= lo) & (Rc < hi)
    if s.sum():
        r = Jm[s] / Jv[s]
        CLASSES[nom] = {"fenetres": int(s.sum()), "J_vrai_median": float(np.median(Jv[s])),
                        "J_modele_median": float(np.median(Jm[s])), "rapport_median": float(np.median(r)),
                        "rapport_min": float(r.min()), "rapport_max": float(r.max())}
        print(f"  {nom:6s} ohms : {s.sum():5d} fenetres ; J_vrai median {np.median(Jv[s]):6.2f}, modele {np.median(Jm[s]):6.2f} ; "
              f"J/J_vrai median {np.median(r):.2f} ({r.min():.2f} a {r.max():.2f})")
cor = float(np.corrcoef(np.log(Jm), np.log(Jv))[0, 1])
print(f"  Correlation de ln J (modele) et ln J_vrai : {cor:.3f} (un jacobien constant donnerait 0).")


# %% ETAPE 6 : export

titre("ETAPE 6 : export")
X_TEST = np.vstack([Xc[i] for i in np.linspace(0, len(Xc) - 1, 5).astype(int)])
Y_TEST, J_TEST = predire(MOD, X_TEST)
savemat(os.path.join(DOSSIER_ELM, "elm_pid_modele.mat"),
        {"W": MOD["W"], "B": MOD["B"].reshape(1, -1), "BETA": MOD["BETA"].reshape(-1, 1), "P0": MOD["P0"],
         "X_MOY": MOD["XM"].reshape(1, -1), "X_EC": MOD["XE"].reshape(1, -1), "T_MOY": MOD["TM"], "T_EC": MOD["TE"],
         "ORDRE": float(ORDRE), "N_CACHES": float(MOD["neurones"]), "C_REGULARISATION": CR["C"],
         "GRAINE": float(CR["premiere_graine"]), "NF": float(NF), "TC": TC, "FORMAT": 2.0,
         "X_TEST": X_TEST, "Y_TEST": Y_TEST.reshape(-1, 1), "J_TEST": J_TEST.reshape(-1, 1),
         "DESCRIPTION": (f"ELM-PID : entrees [ybar(n-1)..ybar(n-{ORDRE}) dbar(n)..dbar(n-{ORDRE})], "
                         f"fenetres de 0.5 ms (110 x Tc), sortie ybar(n) en volts ; J_TEST = d ybar / d dbar(n).")})
with open(os.path.join(DOSSIER_ELM, "entrainement_elm_resultats.json"), "w", encoding="utf-8") as f:
    json.dump({"format": 2, "indice_lipschitz": INDICE, "seuil_coude": SEUIL_COUDE, "ordre": ORDRE,
               "references": REFS, "configurations": {n: {k: v for k, v in c.items()} for n, c in CONFIGS.items()},
               "retenue": RETENUE, "graine": CR["premiere_graine"], "modele": JUGE, "classes_de_charge": CLASSES,
               "correlation_ln_J": cor, "duree_s": round(time.time() - T0, 1)}, f, indent=1, allow_nan=False)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
axes[0].semilogy(list(ORDRES), [INDICE[n] for n in ORDRES], "o-")
axes[0].axvline(ORDRE, color="tab:red", ls="--", lw=0.8)
axes[0].set_xlabel("ordre n")
axes[0].set_title("indice de Lipschitz (He et Asada)", fontsize=9, loc="left")
for nom, lo, hi, col in (("4-6", 4, 6, "tab:blue"), ("6-25", 6, 25, "tab:green"), ("25-50", 25, 50, "tab:orange"),
                         ("50-99", 50, 99, "tab:red")):
    s = (Rc >= lo) & (Rc < hi)
    axes[1].scatter(Jv[s], Jm[s], s=3, color=col, label=f"{nom} ohms")
lim = [0, max(Jv.max(), Jm.max()) * 1.05]
axes[1].plot(lim, lim, "k-", lw=0.6)
axes[1].plot(lim, [x / 2 for x in lim], "k--", lw=0.6)
axes[1].set_xlabel("J vrai (V par unite)")
axes[1].set_ylabel("J du modele")
axes[1].legend(fontsize=7)
axes[1].set_title("jacobien du modele retenu (controle)", fontsize=9, loc="left")
X, T, _ = regresseurs(SIMULINK["ELM_V1"]["y"], SIMULINK["ELM_V1"]["d"], ORDRE)
t = np.arange(len(T))[:300] * 0.5
axes[2].plot(t, T[:300], lw=0.8, label="Simulink")
axes[2].plot(t, predire(MOD, X)[0][:300], lw=0.8, label="prediction a une fenetre")
axes[2].set_xlabel("temps (ms)")
axes[2].legend(fontsize=7)
axes[2].set_title("ELM_V1 : prediction a une fenetre", fontsize=9, loc="left")
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_ELM, "entrainement_elm.png"), dpi=110)
plt.close(fig)
print(f"  elm_pid_modele.mat, entrainement_elm_resultats.json, entrainement_elm.png ecrits ({time.time() - T0:.0f} s).")
