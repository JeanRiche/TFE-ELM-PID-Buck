# =============================================================================
# banc_commun.py
#
# OBJECTIF
#   Banc d'essai Python commun aux cinq methodes. Il simule le meme circuit
#   que le modele Simulink commun Buck_Commun.slx, sur les essais definis
#   par scenarios_communs.json et scenario_*.mat, et calcule les memes
#   grandeurs que les scripts MATLAB Simuler_*.m. Lance seul, il fait tourner
#   le PID de Ziegler-Nichols et ecrit ses resultats attendus : la simulation
#   Simulink de ce PID sur Buck_Commun.slx doit les retrouver. Chaque
#   methode reprend ce banc avec son propre regulateur (la classe du
#   regulateur est le seul element qui change).
#
# LE CIRCUIT SIMULE
#   Buck commute au pas du powergui h = 1/(22000*1200) s (1 200 pas par
#   periode de decoupage), integration par la methode des trapezes, pertes du
#   MOSFET (Ron = 0.1 ohm) et de la diode (Vf = 0.8 V, Ron = 1 mohm). Ce sont
#   les equations du banc qui a predit la boucle ouverte, le PID de
#   Ziegler-Nichols et l'echelon F2 a 0.01 V pres de Simulink. Conventions :
#     - interrupteur : g = 1 si la phase de la porteuse (numero du pas modulo
#       1 200) est sous d x 1 200 ; le pas [j, j+1] utilise l'etat de
#       l'interrupteur aux pas j-1 et j (un pas de retard, comme Simscape) ;
#     - Vin variable (source de tension commandee de Buck_Commun.slx) et
#       charge electronique i = gx x Vout (source de courant commandee), avec
#       la meme convention d'un pas de retard ;
#     - charge = resistance R0 en parallele avec la charge electronique ;
#     - conduction discontinue : interrupteur ouvert et courant de la bobine
#       nul, la diode se bloque ; le courant reste nul et le condensateur se
#       decharge dans la charge.
#   Le regulateur est echantillonne a Tc = 120 pas du powergui. A chaque
#   instant k*Tc : mesure = Vout + bruit, quantifiee au pas q ; consigne =
#   100 + dvref ; rapport cyclique d = sortie du regulateur bornee a
#   [0.01 ; 0.99], applique pendant les 120 pas suivants.
#
# COMMENT LE BANC VA VITE SANS APPROXIMATION
#   Pendant une periode du regulateur, d, Vin et gx sont constants. Les 120
#   pas se regroupent donc en deux suites au plus (MOSFET passant, puis
#   diode passante) ou chaque pas applique la meme application affine a
#   l'etat Z = [iL, Vout, courant de charge du pas precedent, 1]. Une suite
#   de m pas identiques est la puissance m de la matrice d'un pas : le
#   resultat est le meme, au nombre flottant pres, que 120 pas calcules un
#   par un. Les pas ou l'interrupteur ou Vin change sont calcules un par un.
#   En conduction discontinue, l'instant ou le courant s'annule est cherche
#   par dichotomie sur les puissances (le courant decroit de facon monotone
#   quand l'interrupteur est ouvert). Le controle de l'etape 3a compare ce
#   calcul a une simulation pas a pas.
#
# LES GRANDEURS CALCULEES (memes definitions dans les scripts Simuler_*.m)
#   Toutes sont calculees sur les valeurs aux instants k*Tc.
#   e = consigne - Vout (vraie tension, pas la mesure bruitee).
#   Demarrage (de 0 a min(premier evenement, 50 ms)) : temps de montee de
#     10 a 90 V, depassement en % de 100 V, temps d'etablissement dans la
#     bande +-1 V (dernier instant hors bande), erreur moyenne et ondulation
#     crete a crete de 30 ms a la fin du demarrage.
#   Apres le demarrage (de 30 ms a la fin) : IAE, ISE, ITAE (t compte depuis
#     0), ecart maximal, ecart efficace.
#   Pour chaque evenement j (fenetre de t_j a l'evenement suivant) : IAE,
#     ecart maximal, temps de retour dans la bande +-1 V (0 si la tension n'en
#     sort pas) et, sur les 10 dernieres ms de la fenetre, l'amplitude crete
#     a crete de l'ecart (oscillation residuelle) et le fait d'etre revenu
#     dans la bande. Si l'ecart depasse encore 1 V dans ces 10 dernieres ms,
#     la tension n'est "pas revenue" et le temps de retour est laisse vide.
#   Activite de la commande (30 ms a la fin) : moyenne de |d(k) - d(k-1)|,
#     part des instants en butee (d <= 0.01 ou d >= 0.99), ecart type de d.
#   Conduction discontinue : part des instants ou iL < 0.01 A (30 ms a la fin).
#   Courant de la bobine iL (observe seulement) : courant crete pendant le
#     demarrage, ondulation crete a crete de 30 ms a la fin du demarrage,
#     courant maximal et minimal dans la fenetre de chaque evenement. C'est
#     le courant trace sur les figures 6, 13 et 20 de Lu et al. ("output
#     current") : son ondulation y donne L par ses pentes.
#
# FICHIERS
#   Lus : scenario_*.mat et scenarios_communs.json.
#   Ecrits : predictions_banc_pid_classique.json (grandeurs, tension et
#     courant de la bobine toutes les millisecondes, rapport cyclique moyen
#     par milliseconde), banc_commun_pid_classique.png (tension),
#     banc_commun_pid_classique_courant.png (courant de la bobine).
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT :
#   python banc_commun.py
#   Duree : une dizaine de secondes.
#
# =============================================================================


# %% ETAPE 0 : bibliotheques, constantes et lecture des essais

import os                                     # chemins de fichiers
import json                                   # lecture et ecriture au format texte
import time                                   # duree d'execution
import numpy as np                            # calcul numerique
from scipy.io import loadmat                  # lecture des fichiers .mat
import matplotlib.pyplot as plt               # figures

# Circuit (memes valeurs que Buck_Commun.slx)
L_BOB, C_CONV = 10e-3, 47e-6                  # inductance (H), condensateur (F)
RON, VF, RDIODE = 0.1, 0.8, 0.001             # pertes du MOSFET et de la diode
FSW = 22000.0                                 # frequence de decoupage (Hz)
NPER = 1200                                   # pas du powergui par periode de decoupage
H = 1.0 / (FSW * NPER)                        # pas du powergui (s) : 37.88 ns
NREG = 120                                    # pas du powergui par periode du regulateur
TC = NREG * H                                 # periode du regulateur (s) : 1/220000
VREF = 100.0                                  # consigne nominale (V)
D_MIN, D_MAX = 0.01, 0.99                     # saturation physique du rapport cyclique
BANDE = 1.0                                   # bande de retour +-1 V
FIN_FENETRE = 0.010                           # dernieres 10 ms de chaque fenetre d'evenement (s)
SEUIL_DCM = 0.01                              # iL sous 10 mA : conduction discontinue

# PID de Ziegler-Nichols (bloc PID de Buck_Commun.slx) : forme Parallel, discret a Tc,
# Forward Euler pour l'integrateur et le filtre, clamping.
PID_P, PID_I, PID_D, PID_N = 0.093910, 301.089, 7.3227e-06, 64122.9
PID_CI_INTEGRATEUR, PID_CI_FILTRE = 0.5, 0.01

try:                                          # dossier du script
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
print("Dossier de travail :", DOSSIER)

with open(os.path.join(DOSSIER, "scenarios_communs.json"), encoding="utf-8") as f:
    LISTE = json.load(f)                      # liste des essais et choix communs
assert abs(LISTE["Te"] - TC) < 1e-15, "Les essais n'ont pas ete generes a la periode Tc : regenerer les fichiers d'essais a la periode Tc."


def charger_scenario(code):
    """Lit scenario_<code>.mat et renvoie un dictionnaire de tableaux a une dimension."""
    m = loadmat(os.path.join(DOSSIER, f"scenario_{code}.mat"))
    return {"code": code, "nom": str(m["nom"][0]), "duree": float(m["duree"].item()),
            "R0": float(m["R0"].item()), "q": float(m["q"].item()),
            "evenements": [float(x) for x in np.ravel(m["evenements"])],
            "t": m["t"].ravel(), "dvref": m["dvref"].ravel(), "vin": m["vin"].ravel(),
            "gx": m["gx"].ravel(), "bruit": m["bruit"].ravel()}


SCENARIOS = [charger_scenario(e["code"]) for e in LISTE["essais"]]
print(f"{len(SCENARIOS)} essais lus : " + ", ".join(s["code"] for s in SCENARIOS))
# Essais complementaires (S10), hors du cout J. Ils ne sont pas dans SCENARIOS ;
# un script ne les simule que s'il les demande.
SCENARIOS_COMPLEMENTAIRES = [charger_scenario(e["code"]) for e in LISTE.get("essais_complementaires", [])]
if SCENARIOS_COMPLEMENTAIRES:
    print("Essais complementaires lus (hors du cout J) : " + ", ".join(s["code"] for s in SCENARIOS_COMPLEMENTAIRES))


# %% ETAPE 1 : le circuit et la boucle fermee

class Circuit:
    """Un pas du powergui est une application affine de l'etat
    Z = [iL, v, ich_prec, 1], ou ich_prec est le courant de la charge
    electronique au pas precedent. La matrice 4 x 4 d'un pas ne depend que
    de l'interrupteur au debut et a la fin du pas (g0, g1), de gx, de Vin au
    debut et a la fin du pas, et de R0. Les matrices et leurs puissances
    sont gardees en memoire."""

    def __init__(self, R0):
        self.R0 = R0
        self.mat = {}                         # (g0, g1, gx, vin0, vin1) -> matrice d'un pas normal
        self.dcm = {}                         # gx -> matrice d'un pas en conduction discontinue
        self.puis = {}                        # (cle, m) -> puissance m d'une matrice

    def pas_normal(self, g0, g1, gx, vin0, vin1):
        cle = (g0, g1, gx, vin0, vin1)
        if cle not in self.mat:
            a = H / 2.0                                       # demi-pas de la methode des trapezes
            R0 = self.R0
            p0 = -(RON if g0 else RDIODE) / L_BOB             # chute resistive au debut du pas
            p1 = -(RON if g1 else RDIODE) / L_BOB             # et a la fin
            m_fin = np.array([[1.0 - a * p1, a / L_BOB], [-a / C_CONV, 1.0 + a / (R0 * C_CONV)]])
            inv = np.linalg.inv(m_fin)
            # debut du pas ; le courant de charge de fin de pas vaut gx x v(debut)
            m_deb = np.array([[1.0 + a * p0, -a / L_BOB], [a / C_CONV, 1.0 - a / (R0 * C_CONV) - a * gx / C_CONV]])
            s0 = vin0 / L_BOB if g0 else -VF / L_BOB          # source vue par la bobine au debut
            s1 = vin1 / L_BOB if g1 else -VF / L_BOB          # et a la fin du pas
            T = np.zeros((4, 4))
            T[:2, :2] = inv @ m_deb
            T[:2, 2] = inv @ np.array([0.0, -a / C_CONV])     # courant de charge du debut du pas
            T[:2, 3] = inv @ np.array([a * (s0 + s1), 0.0])   # sources
            T[2, 1] = gx                                      # nouveau courant de charge = gx x v
            T[3, 3] = 1.0
            self.mat[cle] = T
        return cle, self.mat[cle]

    def pas_dcm(self, gx):
        if gx not in self.dcm:
            a = H / 2.0
            b = H / (2.0 * self.R0 * C_CONV)                  # terme de decharge du condensateur
            T = np.zeros((4, 4))
            T[1, 1] = (1.0 - b - a * gx / C_CONV) / (1.0 + b)
            T[1, 2] = -a / (C_CONV * (1.0 + b))
            T[2, 1] = gx
            T[3, 3] = 1.0
            self.dcm[gx] = T
        return ("dcm", gx), self.dcm[gx]

    def puissance(self, cle, T, m):
        if m == 0:
            return np.eye(4)
        if (cle, m) not in self.puis:
            self.puis[(cle, m)] = np.linalg.matrix_power(T, m)
        return self.puis[(cle, m)]

    def un_pas(self, Z, g0, g1, gx, vin0, vin1):
        """Un pas, avec la regle de la diode : si l'interrupteur est ouvert
        et que le courant est nul (ou s'annule pendant le pas), le pas est
        calcule en conduction discontinue."""
        if g1 == 0 and Z[0] <= 0.0:
            return self.pas_dcm(gx)[1] @ Z
        Zn = self.pas_normal(g0, g1, gx, vin0, vin1)[1] @ Z
        if g1 == 0 and Zn[0] < 0.0:
            return self.pas_dcm(gx)[1] @ Z
        return Zn

    def suite(self, Z, g, n, gx, vin):
        """n pas identiques (interrupteur g au debut et a la fin de chaque pas)."""
        if n == 0:
            return Z
        cle, T = self.pas_normal(g, g, gx, vin, vin)
        if g == 1:
            return self.puissance(cle, T, n) @ Z
        cd, Td = self.pas_dcm(gx)
        if Z[0] <= 0.0:                                       # deja en conduction discontinue
            return self.puissance(cd, Td, n) @ Z
        Zn = self.puissance(cle, T, n) @ Z
        if Zn[0] >= 0.0:                                      # le courant ne s'annule pas
            return Zn
        bas, haut = 1, n                                      # premier pas ou le courant passerait sous 0
        while bas < haut:
            mil = (bas + haut) // 2
            if (self.puissance(cle, T, mil) @ Z)[0] < 0.0:
                haut = mil
            else:
                bas = mil + 1
        Zm = self.puissance(cle, T, bas - 1) @ Z              # etat avant ce pas (courant encore positif)
        Zm = Td @ Zm                                          # ce pas est calcule en conduction discontinue
        return self.puissance(cd, Td, n - bas) @ Zm


def simuler(regulateur, sc):
    """Boucle fermee sur l'essai sc. A chaque instant k*Tc :
      1. lecture de Vout ; mesure = Vout + bruit, quantifiee au pas q ;
      2. erreur = (100 + dvref) - mesure ; commande mu du regulateur ;
      3. rapport cyclique d = mu borne a [0.01, 0.99] ;
      4. 120 pas du circuit avec ce rapport cyclique, Vin(k) et gx(k).
    Un regulateur qui a besoin de la mesure (attribut utilise_mesure) la
    recoit en second argument : c'est la sortie du capteur, jamais 100 - e."""
    circ = Circuit(sc["R0"])
    N = len(sc["t"])
    sortie = {k: np.zeros(N) for k in ("v", "mesure", "consigne", "mu", "d", "iL")}
    Z = np.array([0.0, 0.0, 0.0, 1.0])                        # circuit au repos
    g_prec, vin_prec = 0, float(sc["vin"][0])
    for k in range(N):
        v = Z[1]
        mesure = v + sc["bruit"][k]
        mesure = sc["q"] * np.round(mesure / sc["q"])          # quantification (pas de 1e-9 V : sans effet)
        consigne = VREF + sc["dvref"][k]
        e = consigne - mesure
        mu = regulateur.pas(e, mesure) if getattr(regulateur, "utilise_mesure", False) else regulateur.pas(e)
        d = min(max(mu, D_MIN), D_MAX)
        sortie["v"][k], sortie["mesure"][k], sortie["consigne"][k] = v, mesure, consigne
        sortie["mu"][k], sortie["d"][k], sortie["iL"][k] = mu, d, Z[0]
        if k == N - 1:
            break
        gx, vin = float(sc["gx"][k]), float(sc["vin"][k])
        ph0 = (k * NREG) % NPER                               # phase de la porteuse au debut de la periode
        n_on = int(np.ceil(d * NPER))                         # pas de phase ph < d x 1200 dans une periode
        n1 = min(max(n_on - ph0, 0), NREG)                    # pas MOSFET passant dans cette periode
        for g, n in ((1, n1), (0, NREG - n1)):
            if n == 0:
                continue
            if g != g_prec or vin != vin_prec:                # premier pas calcule a part
                Z = circ.un_pas(Z, g_prec, g, gx, vin_prec, vin)
                n -= 1
                g_prec, vin_prec = g, vin
            Z = circ.suite(Z, g, n, gx, vin)
    sortie["t"] = sc["t"]
    return sortie


# %% ETAPE 2 : les grandeurs mesurees (memes definitions que les scripts Simuler_*.m)

def grandeurs(sim, sc):
    """Toutes les grandeurs d'un essai, en numeros d'instants k (pas Tc)."""
    v, consigne, d, iL = sim["v"], sim["consigne"], sim["d"], sim["iL"]
    N = len(v)
    n = np.arange(N)
    e = consigne - v
    k30 = int(round(0.03 / TC))
    evts = [int(round(x / TC)) for x in sc["evenements"]]
    k_dem = min([int(round(0.05 / TC))] + evts)                 # fin du demarrage
    o = {}
    # Demarrage
    k10 = np.where(v >= 10.0)[0]
    k90 = np.where(v >= 90.0)[0]
    o["t_montee_ms"] = float((k90[0] - k10[0]) * TC * 1e3) if len(k10) and len(k90) else float("nan")
    o["depassement_pct"] = float(max(0.0, (v[:k_dem].max() - VREF) / VREF * 100.0))
    hors = np.where(np.abs(e[:k_dem]) > BANDE)[0]
    o["t_etab_ms"] = float((hors[-1] + 1) * TC * 1e3) if len(hors) else 0.0
    o["erreur_moy_dem_V"] = float(np.mean(e[k30:k_dem]))
    o["ondulation_V"] = float(v[k30:k_dem].max() - v[k30:k_dem].min())
    o["iL_max_dem_A"] = float(iL[:k_dem].max())                 # courant crete du demarrage
    o["iL_ondulation_A"] = float(iL[k30:k_dem].max() - iL[k30:k_dem].min())   # ondulation en regime
    # Apres le demarrage
    w = n >= k30
    o["IAE"] = float(np.sum(np.abs(e[w])) * TC)
    o["ISE"] = float(np.sum(e[w] ** 2) * TC)
    o["ITAE"] = float(np.sum(n[w] * TC * np.abs(e[w])) * TC)
    o["e_max_V"] = float(np.max(np.abs(e[w])))
    o["e_eff_V"] = float(np.sqrt(np.mean(e[w] ** 2)))
    # Evenements
    o["evenements"] = []
    bornes = evts + [N]
    n_fin = int(round(FIN_FENETRE / TC))                          # 2 200 instants
    for j, kj in enumerate(evts):
        fen = (n >= kj) & (n < bornes[j + 1])                     # de l'evenement au suivant
        fin = (n >= max(kj, bornes[j + 1] - n_fin)) & (n < bornes[j + 1])   # ses 10 dernieres ms
        hors = np.where(fen & (np.abs(e) > BANDE))[0]             # instants hors de la bande +-1 V
        revenu = bool(np.max(np.abs(e[fin])) <= BANDE)            # dans la bande sur les 10 dernieres ms
        if len(hors) == 0:
            t_ret = 0.0                                           # jamais sorti de la bande
        elif revenu:
            t_ret = float((hors[-1] - kj + 1) * TC * 1e3)         # dernier instant hors bande, plus un pas
        else:
            t_ret = None                                          # pas revenu avant l'evenement suivant
        o["evenements"].append({"t_ms": kj * TC * 1e3,
                                "IAE": float(np.sum(np.abs(e[fen])) * TC),
                                "e_max_V": float(np.max(np.abs(e[fen]))),
                                "t_retour_ms": t_ret,
                                "reste_dans_bande": bool(len(hors) == 0),
                                "revenu": revenu,
                                "e_crete_a_crete_fin_V": float(np.ptp(e[fin])),
                                "iL_max_A": float(np.max(iL[fen])),      # courant observe dans la fenetre
                                "iL_min_A": float(np.min(iL[fen]))})
    # Activite de la commande et conduction discontinue
    dw = d[w]
    o["dd_moyen"] = float(np.mean(np.abs(np.diff(dw))))
    o["butee_pct"] = float(np.mean((dw <= D_MIN + 1e-12) | (dw >= D_MAX - 1e-12)) * 100.0)
    o["d_ecart_type"] = float(np.std(dw))
    o["dcm_pct"] = float(np.mean(iL[w] < SEUIL_DCM) * 100.0)
    return o


def resume(o):
    """Une ligne lisible des grandeurs principales."""
    return (f"montee {o['t_montee_ms']:.2f} ms, depassement {o['depassement_pct']:.2f} %, etabli a "
            f"{o['t_etab_ms']:.2f} ms ; IAE {o['IAE']*1e3:.2f} mV.s, ecart max {o['e_max_V']:.2f} V, "
            f"efficace {o['e_eff_V']*1e3:.1f} mV ; |dd| {o['dd_moyen']:.4f}, butee {o['butee_pct']:.1f} %, "
            f"DCM {o['dcm_pct']:.1f} % ; iL crete au demarrage {o['iL_max_dem_A']:.1f} A")


# %% ETAPE 3 : le PID classique, et les controles du banc

class PIDClassique:
    """Copie du bloc "PID Controller" de Buck_Commun.slx : forme
    Parallel, Forward Euler a Tc pour l'integrateur et le filtre, sortie
    dans [0.01, 0.99], clamping : l'integration s'arrete quand la somme des
    termes sort des bornes et que l'entree de l'integrateur a le meme signe
    que le depassement."""

    def __init__(self):
        self.xI = PID_CI_INTEGRATEUR          # etat de l'integrateur
        self.xF = PID_CI_FILTRE               # etat du filtre de derivee (D x e filtre)

    def pas(self, e):
        derivee = PID_N * (PID_D * e - self.xF)          # terme derive filtre
        u = PID_P * e + self.xI + derivee                # somme des trois termes
        u_sat = min(max(u, D_MIN), D_MAX)
        entree_I = PID_I * e
        bloque = (u != u_sat) and (np.sign(entree_I) == np.sign(u - u_sat))
        if not bloque:
            self.xI += TC * entree_I
        self.xF += TC * derivee
        return u_sat


def simuler_pas_a_pas(regulateur, sc, duree):
    """Reference lente : le meme circuit calcule pas du powergui par pas,
    sans regroupement. Sert seulement au controle de l'etape 3a."""
    circ = Circuit(sc["R0"])
    K = int(round(duree / TC))
    Z = np.array([0.0, 0.0, 0.0, 1.0])
    g_prec, vin_prec = 0, float(sc["vin"][0])
    vs = np.zeros(K)
    for k in range(K):
        v = Z[1]
        mesure = sc["q"] * np.round((v + sc["bruit"][k]) / sc["q"])
        d = min(max(regulateur.pas(VREF + sc["dvref"][k] - mesure), D_MIN), D_MAX)
        vs[k] = v
        gx, vin = float(sc["gx"][k]), float(sc["vin"][k])
        for m in range(NREG):
            j = k * NREG + m
            g = 1 if (j % NPER) < d * NPER else 0
            Z = circ.un_pas(Z, g_prec, g, gx, vin_prec, vin)
            g_prec, vin_prec = g, vin
    return vs


# === FIN DES DEFINITIONS DU BANC COMMUN ===
# Les scripts qui reprennent ce banc (recherche_pid_fige.py, puis ceux des
# methodes) executent ce fichier jusqu'a la ligne ci-dessus : constantes,
# lecture des essais, circuit, boucle fermee, grandeurs et PID classique.
# Ce qui suit ne tourne que lorsqu'on lance banc_commun.py lui-meme.

T0 = time.time()
titre = lambda txt: print("\n" + "=" * 78 + "\n " + txt + "\n" + "=" * 78)

titre("ETAPE 3a : controles du banc")
S1 = SCENARIOS[0]
# 1. Regroupement des pas contre calcul pas a pas, sur un essai qui fait tout
#    varier : charge electronique a 98 ohms + conductance (5 ohms), echelon de
#    Vin a 2 ms, delestage vers 98 ohms a 3 ms (conduction discontinue).
n_test = int(round(0.006 / TC)) + 1
test = dict(S1, R0=98.0, t=S1["t"][:n_test], dvref=np.zeros(n_test), bruit=np.zeros(n_test),
            vin=np.where(np.arange(n_test) * TC < 0.002, 200.0, 160.0),
            gx=np.where(np.arange(n_test) * TC < 0.003, 1.0 / 5.0 - 1.0 / 98.0, 0.0))
v_rapide = simuler(PIDClassique(), test)["v"][:n_test - 1]
v_lent = simuler_pas_a_pas(PIDClassique(), test, (n_test - 1) * TC)
print(f"  Regroupement des pas contre calcul pas a pas (6 ms, Vin et charge variables, conduction "
      f"discontinue) : ecart max {np.max(np.abs(v_rapide - v_lent)):.2e} V.")
# 2. Grandeurs de reference mesurees sous Simulink sur Buck_Commun.slx
sim_ref = simuler(PIDClassique(), S1)
tm = sim_ref["v"][: (len(sim_ref["v"]) // 10) * 10].reshape(-1, 10).mean(axis=1)   # moyenne par periode
t99 = np.argmax(tm >= 99.0) * 10 * TC * 1e3
print(f"  S1 : 99 V (moyenne par periode) a {t99:.3f} ms, maximum moyen {tm[: int(0.05 / TC / 10)].max():.2f} V "
      f"(Simulink, Buck_Commun.slx : 1.545 ms et 113.96 V).")
# 3. Charge electronique : 98 ohms en parallele avec gx = 1/5 - 1/98 doit valoir 5 ohms
S1_elec = dict(S1, R0=98.0, gx=np.full(len(S1["t"]), 1.0 / 5.0 - 1.0 / 98.0))
sim_elec = simuler(PIDClassique(), S1_elec)
ecart = np.abs(sim_elec["v"] - sim_ref["v"])
print(f"  Charge electronique contre resistance de 5 ohms : ecart max {ecart.max()*1e3:.1f} mV, "
      f"efficace {np.sqrt(np.mean(ecart**2))*1e3:.1f} mV.")
controle = {"regroupement_ecart_max_V": float(np.max(np.abs(v_rapide - v_lent))),
            "S1_99V_moyenne_ms": float(t99), "S1_max_moyen_V": float(tm[: int(0.05 / TC / 10)].max()),
            "S1_99V_ms": float(np.argmax(sim_ref["v"] >= 99.0) * TC * 1e3), "S1_max_V": float(sim_ref["v"].max()),
            "charge_electronique_ecart_max_V": float(ecart.max())}

titre("ETAPE 3b : le PID classique sur les onze essais (puis les essais complementaires, hors du cout J)")
predictions = {"format": 2, "Te": TC, "regulateur": "PID de Ziegler-Nichols (Buck_Commun.slx)",
               "controle": controle, "essais": {}}
sims = {}
for sc in SCENARIOS + SCENARIOS_COMPLEMENTAIRES:
    t0 = time.time()
    sim = simuler(PIDClassique(), sc)
    o = grandeurs(sim, sc)
    sims[sc["code"]] = sim
    pas_1ms = int(round(1e-3 / TC))
    predictions["essais"][sc["code"]] = {"grandeurs": o,
                                         "v_toutes_les_ms": sim["v"][::pas_1ms].tolist(),
                                         "iL_toutes_les_ms": sim["iL"][::pas_1ms].tolist(),
                                         "d_moyen_par_ms": [float(np.mean(sim["d"][i:i + pas_1ms]))
                                                            for i in range(0, len(sim["d"]) - pas_1ms + 1, pas_1ms)]}
    print(f"  {sc['code']:4s} {sc['nom']:24s} {resume(o)}  [{time.time()-t0:.0f} s]")
    for ev in o["evenements"]:
        if ev["reste_dans_bande"]:
            etat = "reste dans la bande"
        elif ev["revenu"]:
            etat = f"retour en {ev['t_retour_ms']:.2f} ms"
        else:
            etat = "PAS REVENU dans la bande avant l'evenement suivant"
        print(f"         evenement a {ev['t_ms']:6.1f} ms : IAE {ev['IAE']*1e3:7.2f} mV.s, ecart max {ev['e_max_V']:6.2f} V, "
              f"ecart crete a crete en fin de fenetre {ev['e_crete_a_crete_fin_V']:5.2f} V, {etat} ; "
              f"iL de {ev['iL_min_A']:.1f} a {ev['iL_max_A']:.1f} A")


# %% ETAPE 4 : fichier des resultats attendus et figure

with open(os.path.join(DOSSIER, "predictions_banc_pid_classique.json"), "w", encoding="utf-8") as f:
    json.dump(predictions, f, indent=1, allow_nan=False)
print("\n  predictions_banc_pid_classique.json ecrit.")

fig, axes = plt.subplots(len(SCENARIOS), 1, figsize=(11, 2.0 * len(SCENARIOS)))
for ax, sc in zip(axes, SCENARIOS):
    sim = sims[sc["code"]]
    ax.plot(sim["t"] * 1e3, sim["v"], lw=0.6, label="Vout")
    ax.plot(sim["t"] * 1e3, sim["consigne"], "k--", lw=0.6, label="consigne")
    ax.set_ylim(min(40.0, sim["consigne"].min() - 10), max(125.0, sim["consigne"].max() + 15))
    ax.set_title(f"{sc['code']} : {sc['nom']} (PID classique)", fontsize=9, loc="left")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, loc="lower right")
axes[-1].set_xlabel("temps (ms)")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER, "banc_commun_pid_classique.png"), dpi=110)
plt.close(fig)
print("  banc_commun_pid_classique.png ecrit.")

# Courant de la bobine sur les onze essais (observe seulement)
fig, axes = plt.subplots(len(SCENARIOS), 1, figsize=(11, 2.0 * len(SCENARIOS)))
for ax, sc in zip(axes, SCENARIOS):
    sim = sims[sc["code"]]
    ax.plot(sim["t"] * 1e3, sim["iL"], lw=0.5, color="tab:red", label="iL")
    ax.set_ylabel("A")
    ax.set_title(f"{sc['code']} : {sc['nom']}, courant de la bobine (PID classique)", fontsize=9, loc="left")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, loc="upper right")
axes[-1].set_xlabel("temps (ms)")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER, "banc_commun_pid_classique_courant.png"), dpi=110)
plt.close(fig)
print(f"  banc_commun_pid_classique_courant.png ecrit. Duree totale : {time.time() - T0:.0f} s.")
