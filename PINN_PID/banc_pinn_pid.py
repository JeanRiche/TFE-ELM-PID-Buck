# =============================================================================
# banc_pinn_pid.py
#
# OBJECTIF
#   Faire tourner le PINN-PID (correction retenue par choix_correction_pinn.py,
#   lue dans choix_correction_pinn.json) sur le banc commun, sur les onze
#   essais de developpement et sur S10, a cote des references (PID de
#   Ziegler-Nichols, qui est aussi le PINN-PID sans adaptation, et meilleur
#   PID fige sur grille), et mesurer ce que quelques pieces apportent
#   (ablations).
#
# LE PINN-PID, TEL QU'IL EST CODE ICI
#   A chaque periode Tc = 1/220 000 s, le bloc PID de la base commune
#   calcule u avec les gains K = K_ZN .* x ; observateur (filtre de Kalman
#   etendu) ; a la fin de chaque fenetre de 110 periodes (0.5 ms),
#   5 iterations d'Adam (reglages d'Ito et Wasa) sur x, depart a chaud,
#   moments remis a zero, horizon de 110 pas avec le bloc PID et le PINN
#   dans la boucle, gradient par retropropagation dans le temps.
#   Cout J = moyenne de 1/2 (e^2 + W_U u^2) + W_F |x - c|^2.
#   Corrections candidates :
#     C1  zone morte : l'optimisation n'a lieu que si l'erreur EFFICACE de
#         la fenetre depasse 0.1 V (Peterson et Narendra 1982 ; Ioannou et
#         Sun 1996) ; c'est la correction retenue ;
#     C2  regularisation vers les gains de Ziegler-Nichols : c = (1, 1, 1)
#         au lieu de c = 0 (Tikhonov vers un a priori), W_F inchange ;
#     C3  ensemble admissible de l'ELM-PID (ensemble_gains_elm.mat, copie
#         identique de ELM_PID) et projection de l'ELM-PID (fonctions
#         g_interp, gradient_g, restaurer, projeter recopiees a l'identique
#         de ELM_PID/banc_elm_pid.py) au lieu de la boite.
#   Sans correction (ni zone morte), J = 0.893 sur les onze essais
#   (ablation).
#
# CONTROLES (avant tout calcul sur les onze essais)
#   - gradient de l'horizon contre differences finies (c = 0 et c = 1) ;
#   - fonctions de projection identiques au texte de ELM_PID/banc_elm_pid.py
#     et table identique octet pour octet (si le dossier ELM_PID est la) ;
#   - sans adaptation, ecart 0 sur d au PID de Ziegler-Nichols (S1, S8b).
#
# FICHIERS NECESSAIRES (dans le meme dossier que ce script)
#   banc_commun.py, scenarios_communs.json, scenario_S*.mat,
#   predictions_banc_pid_fige.json, boite_gains_pinn.json,
#   ensemble_gains_elm.mat, pinn_pid_modele.mat, choix_correction_pinn.json
#   (ecrit par choix_correction_pinn.py).
#
# CE QUE PRODUIT CE SCRIPT
#   banc_pinn_pid_resultats.json, banc_pinn_pid.png, et ce que Simulink doit
#   retrouver (meme format que l'ELM-PID, ELM_PID/banc_elm_pid.py) :
#     - predictions_banc_pinn_pid.json (format 2, lu par Simuler_PINN_PID.m
#       et les scripts de construction) ;
#     - pinn_pid_reglages.mat (reglages lus par pinn_pid_adaptatif.m,
#       FORMAT = 3) ;
#     - reference_rejeu_pinn_pid.mat (lu par Tester_PINN_PID_Rejeu.m : e,
#       mesure, u, gains pas par pas et decisions par fenetre sur S1, S3,
#       S7b, S8b, S9).
#   La classe AdaptateurPINN est la copie Python du bloc MATLAB
#   pinn_pid_adaptatif.m : memes entrees (e, mesure, u), meme sortie (gains),
#   meme ordre des operations. Elle refait elle-meme les etats du bloc PID
#   (integrateur et filtre) a partir de e et de ses propres gains ; PINNPID
#   verifie a chaque pas que cette copie est identique au bloc.
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT
#   python banc_pinn_pid.py      (quatre processus ; dix a vingt minutes)
#
# ORDRE D'EXECUTION
#   1. choix_correction_pinn.py (controles, essais de reglage E1 a E4, choix
#   mecanique de la correction) ; 2. ce script ; puis, dans MATLAB :
#   3. Tester_PINN_PID_Rejeu.m ; 4. Construction_PINN_PID.m ;
#   5. Simuler_PINN_PID.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques, reglages et definitions du banc commun

import os                                     # chemins de fichiers
import json                                   # lecture et ecriture au format texte
import time                                   # duree d'execution
import types                                  # resultats legers renvoyes par les processus
import multiprocessing as mp                  # essais en parallele (un processus par essai)
import numpy as np                            # calcul numerique
from scipy.io import loadmat, savemat         # fichiers .mat
from scipy.linalg import expm                 # table du modele physique
from scipy.signal import cont2discrete        # marges du modele moyen
import matplotlib
matplotlib.use("Agg")                         # figure ecrite dans un fichier, sans fenetre
import matplotlib.pyplot as plt               # figure

try:                                          # dossier du script (ou dossier courant)
    DOSSIER_PINN = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PINN = os.getcwd()

with open(os.path.join(DOSSIER_PINN, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC or "iL_max_dem_A" not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py ne contient pas les definitions attendues (marque de fin, observation du courant).")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])  # constantes, essais, Circuit, simuler, grandeurs, PIDClassique

T_DEBUT = time.time()                         # pour la duree totale
ESSAIS = {sc["code"]: sc for sc in SCENARIOS}  # les onze essais, par code
CODES_IAE = ["S1", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]   # IAE de 30 ms a la fin
CODES_DEM = ["S1", "S8a", "S8b"]              # IAE du demarrage (5, 25 et 98 ohms)
CODES_REJEU = ["S1", "S3", "S7b", "S8b", "S9"]   # essais du rejeu pas a pas (Tester_PINN_PID_Rejeu.m)
K30 = int(round(0.03 / TC))                   # instant 30 ms (numero de pas Tc)
NF = 110                                      # periodes Tc par fenetre de 0.5 ms
TCN = TC * PID_N                              # Tc x N du filtre de derivee
L_MOD, C_MOD = L_BOB, C_CONV                  # L et C du modele du regulateur (le circuit peut etre perturbe a part)
K_ZN = np.array([PID_P, PID_I, PID_D])       # Ziegler-Nichols, gains du bloc PID [P, I, D]

# Boite des gains (boite_gains_pinn.py : marge et coupure aux coins nominaux, loi du bloc parallele)
with open(os.path.join(DOSSIER_PINN, "boite_gains_pinn.json"), encoding="utf-8") as f:
    BOITE = json.load(f)
if np.max(np.abs(np.array(BOITE["K_depart"]) - K_ZN) / K_ZN) > 1e-12:
    raise RuntimeError("boite_gains_pinn.json n'a pas ete calcule avec les gains de Ziegler-Nichols de ce banc.")
X_MIN = np.array(BOITE["multiplicateurs_min"])   # bornes de la boite en multiplicateurs de ZN
X_MAX = np.array(BOITE["multiplicateurs_max"])

# Ensemble admissible de l'ELM-PID (correction C3) : table identique a ELM_PID/ensemble_gains_elm.mat
ENS = loadmat(os.path.join(DOSSIER_PINN, "ensemble_gains_elm.mat"))
G_TABLE = ENS["G_TABLE"].astype(float)
LOG2_MIN, PAS_LOG2 = float(ENS["LOG2_MIN"].item()), float(ENS["PAS_LOG2"].item())
N_TABLE = G_TABLE.shape[0]
if np.max(np.abs(ENS["K_DEPART"].ravel() - K_ZN) / K_ZN) > 1e-12:
    raise RuntimeError("ensemble_gains_elm.mat n'a pas ete calcule avec les gains de Ziegler-Nichols de ce banc.")
# Reglages de la projection : ceux de l'ELM-PID (ELM_PID/banc_elm_pid.py, REGLAGES)
R_PROJ = {"MULT_MIN": 0.25, "MULT_MAX": 4.0, "N_BISSECTIONS": 30.0, "H_GRADIENT": 1e-4,
          "GLISSEMENT": 1.0, "RESTAURATION": 1.0, "N_PROJECTION": 10.0}

# >>> DEBUT DE LA COPIE de ELM_PID/banc_elm_pid.py (de "def g_interp" a "class AdaptateurELM", sans changement)
def g_interp(x):
    """g de l'ensemble admissible, interpole (trilineaire en log2 des
    multiplicateurs). Memes operations que g_interp du bloc MATLAB."""
    u = (np.log2(x) - LOG2_MIN) / PAS_LOG2
    u = np.minimum(np.maximum(u, 0.0), N_TABLE - 1.0)
    i = np.minimum(np.floor(u), N_TABLE - 2.0).astype(int)
    f = u - i
    s = 0.0
    for di in (0, 1):
        wi = f[0] if di else 1.0 - f[0]
        for dj in (0, 1):
            wj = f[1] if dj else 1.0 - f[1]
            for dk in (0, 1):
                wk = f[2] if dk else 1.0 - f[2]
                s += wi * wj * wk * G_TABLE[i[0] + di, i[1] + dj, i[2] + dk]
    return s


def gradient_g(x, h):
    """Gradient de g par differences centrees (normale a la frontiere)."""
    gr = np.zeros(3)
    for q in range(3):
        ep = np.zeros(3)
        ep[q] = h
        gr[q] = (g_interp(x + ep) - g_interp(x - ep)) / (2.0 * h)
    return gr


def restaurer(p, gr, r):
    """Restauration : depuis p (hors de l'ensemble), recul le long de
    -gr jusqu'a la frontiere (doublement puis bissection) ; None si un recul
    de 1 ne suffit pas."""
    lo_, hi_ = r["MULT_MIN"], r["MULT_MAX"]
    nh = gr / np.sqrt(gr @ gr)
    s_ok = 1e-6
    while g_interp(np.minimum(np.maximum(p - s_ok * nh, lo_), hi_)) > 0.0:
        s_ok *= 2.0
        if s_ok > 1.0:
            return None
    s_ko = 0.0 if s_ok == 1e-6 else 0.5 * s_ok
    for _ in range(int(r["N_BISSECTIONS"])):
        m = 0.5 * (s_ko + s_ok)
        if g_interp(np.minimum(np.maximum(p - m * nh, lo_), hi_)) <= 0.0:
            s_ok = m
        else:
            s_ko = m
    return np.minimum(np.maximum(p - s_ok * nh, lo_), hi_)


def projeter(xc, dx, r):
    """Pas dx depuis xc (admissible), projete sur l'ensemble admissible.
    RESTAURATION = 1 : point admissible le plus proche de la cible
    z = xc + dx, par la methode de projection du gradient de Rosen (1961) :
    depuis le point de sortie, on retire la composante sortante de (z - y)
    le long de la normale, on avance dans le plan tangent, puis on revient
    sur la frontiere le long de la normale (restauration) ; on recommence
    tant que le point se rapproche de z (N_PROJECTION fois au plus).
    RESTAURATION = 0 : glissement puis bissection le long du pas tangent
    (variante conservee pour l'ablation).
    Sans glissement : arret au point de sortie."""
    lo_, hi_ = r["MULT_MIN"], r["MULT_MAX"]
    nb = int(r["N_BISSECTIONS"])
    xn = np.minimum(np.maximum(xc + dx, lo_), hi_)
    if g_interp(xn) <= 0.0:
        return xn
    z = xn
    a, b = 0.0, 1.0
    for _ in range(nb):
        m = 0.5 * (a + b)
        if g_interp(np.minimum(np.maximum(xc + m * dx, lo_), hi_)) <= 0.0:
            a = m
        else:
            b = m
    xb = np.minimum(np.maximum(xc + a * dx, lo_), hi_)
    if not r["GLISSEMENT"]:
        return xb
    if r["RESTAURATION"]:
        y = xb
        for _ in range(int(r["N_PROJECTION"])):
            gr = gradient_g(y, r["H_GRADIENT"])
            if gr @ gr == 0.0:
                break
            d = z - y
            sortant = gr @ d
            if sortant > 0.0:
                d = d - sortant / (gr @ gr) * gr
            pt = np.minimum(np.maximum(y + d, lo_), hi_)
            if g_interp(pt) > 0.0:
                pt = restaurer(pt, gr, r)
                if pt is None:
                    break
            if np.sum((z - pt) ** 2) >= np.sum((z - y) ** 2):
                break
            y = pt
        return y
    gr = gradient_g(xb, r["H_GRADIENT"])
    reste = (1.0 - a) * dx
    sortant = gr @ reste
    if sortant > 0.0:
        reste = reste - sortant / (gr @ gr) * gr
    xn = np.minimum(np.maximum(xb + reste, lo_), hi_)
    if g_interp(xn) <= 0.0:
        return xn
    a, b = 0.0, 1.0
    for _ in range(nb):
        m = 0.5 * (a + b)
        if g_interp(np.minimum(np.maximum(xb + m * reste, lo_), hi_)) <= 0.0:
            a = m
        else:
            b = m
    return np.minimum(np.maximum(xb + a * reste, lo_), hi_)
# <<< FIN DE LA COPIE


def _bloc_copie(texte):
    a = texte.index("# >>> DEBUT DE LA COPIE")
    a = texte.index("\n", a) + 1
    return texte[a:texte.index("# <<< FIN DE LA COPIE")]


# Controle C3 : meme table (octet pour octet) et memes fonctions (texte) que l'ELM-PID, si son dossier est la
_DOSSIER_ELM = os.path.join(os.path.dirname(DOSSIER_PINN), "ELM_PID")
if os.path.isdir(_DOSSIER_ELM):
    with open(os.path.join(_DOSSIER_ELM, "ensemble_gains_elm.mat"), "rb") as f1, \
            open(os.path.join(DOSSIER_PINN, "ensemble_gains_elm.mat"), "rb") as f2:
        if f1.read() != f2.read():
            raise RuntimeError("ensemble_gains_elm.mat differe de celui de ELM_PID.")
    with open(os.path.join(_DOSSIER_ELM, "banc_elm_pid.py"), encoding="utf-8") as f1:
        _t = f1.read()
    with open(os.path.join(DOSSIER_PINN, "banc_pinn_pid.py"), encoding="utf-8") as f2:
        _copie = _bloc_copie(f2.read())
    if _t[_t.index("def g_interp(x):"):_t.index("class AdaptateurELM:")].rstrip() + "\n" != _copie:
        raise RuntimeError("Les fonctions de projection different de celles de ELM_PID/banc_elm_pid.py.")
    _ELM_IDENTIQUE = "table et fonctions identiques a ELM_PID"
else:
    _ELM_IDENTIQUE = "dossier ELM_PID absent : comparaison non faite"
if os.path.isdir(_DOSSIER_ELM):                # memes reglages de projection que l'ELM-PID (valeurs lues)
    import re
    for _k, _v in R_PROJ.items():
        _m = re.search(r'"' + _k + r'": ([^,}]+)', _t)
        if _m is None or float(_m.group(1)) != _v:
            raise RuntimeError(f"Reglage de projection {_k} different de celui de l'ELM-PID.")
if g_interp(np.ones(3)) > 0.0:
    raise RuntimeError("Le depart (gains de Ziegler-Nichols) n'est pas admissible dans la table.")
print(f"Ensemble admissible de l'ELM-PID (C3) : {N_TABLE}^3 points, {np.mean(G_TABLE <= 0):.1%} admissibles ; "
      f"{_ELM_IDENTIQUE}.")

# Reglages du PINN-PID ; SEUIL_C1 : zone morte de la correction C1
REGLAGES = {"NH": float(NF), "ALPHA": 0.01, "BETA1": 0.9, "BETA2": 0.999, "EPS_ADAM": 1e-7, "N_ADAM": 5.0,
            "W_U": 1e-5, "W_F": 1e-3, "PENTE_SAT": 0.01, "SEUIL_C1": 0.1,
            "G_NOM": 0.2, "VIN_DEPART": 200.0, "SIG_Y": (200.0 / 4096.0) / np.sqrt(12.0),
            "SIG_V": 1e-3, "SIG_I": 1e-2, "SIG_VIN": 0.03, "P0_I": 1.0, "P0_VIN": 30.0 ** 2}

FICHIER_CHOIX = os.path.join(DOSSIER_PINN, "choix_correction_pinn.json")


def combinaison_retenue():
    """Combinaison de corrections choisie mecaniquement sur E1 a E4 par
    choix_correction_pinn.py."""
    if not os.path.exists(FICHIER_CHOIX):
        raise RuntimeError("choix_correction_pinn.json absent : lancer d'abord choix_correction_pinn.py.")
    with open(FICHIER_CHOIX, encoding="utf-8") as f:
        return tuple(json.load(f)["combinaison_retenue"])


# Meilleur PID fige (deuxieme reference, oracle)
with open(os.path.join(DOSSIER_PINN, "predictions_banc_pid_fige.json"), encoding="utf-8") as f:
    GAINS_FIGE = json.load(f)["gains"]

# Modele PINN (entrainement_pinn.py)
MODELE = loadmat(os.path.join(DOSSIER_PINN, "pinn_pid_modele.mat"))
COUCHES = []
l = 1
while f"W{l}" in MODELE:
    COUCHES.append((MODELE[f"W{l}"].astype(float), MODELE[f"b{l}"].astype(float).ravel()))
    l += 1
if abs(float(MODELE["TC"].item()) - TC) > 1e-18:
    raise RuntimeError("pinn_pid_modele.mat n'a pas ete appris avec la periode Tc de ce banc.")
SV_P, SI_P = float(MODELE["SV"].item()), float(MODELE["SI"].item())          # echelles de sortie du PINN
G_C, G_E = float(MODELE["G_CENTRE"].item()), float(MODELE["G_ECHELLE"].item())   # normalisation de G
print(f"Boite (multiplicateurs de ZN) : P de {X_MIN[0]:.3f} a {X_MAX[0]:.3f}, I de {X_MIN[1]:.3f} a {X_MAX[1]:.3f}, "
      f"D de {X_MIN[2]:.3f} a {X_MAX[2]:.3f} ; depart a Ziegler-Nichols (1, 1, 1).")
print(f"PINN : {len(COUCHES) - 1} couches cachees de {COUCHES[0][0].shape[1]} neurones "
      f"(candidat {str(MODELE['CANDIDAT'][0])}, plage {str(MODELE['PLAGE'][0])}).")


def titre(texte):
    """Bandeau de titre dans la console."""
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78, flush=True)


# %% ETAPE 1 : predicteurs, observateur, horizon et regulateur

# Table du modele physique moyen par tranche a charge nominale (inchangee) :
#   x(fin) = PHI(s) x(debut) + GAM(s) (s Vin - (1 - s) Vf) / L.
NS = NREG
PHI = np.zeros((NS + 1, 2, 2)); GAM = np.zeros((NS + 1, 2))
for j in range(NS + 1):
    s = j / NS
    req = s * RON + (1 - s) * RDIODE          # resistance moyenne vue par la bobine
    A = np.array([[-REGLAGES["G_NOM"] / C_MOD, 1 / C_MOD], [-1 / L_MOD, -req / L_MOD]])
    M = np.zeros((4, 4)); M[:2, :2] = A * TC; M[:2, 2:] = np.eye(2) * TC
    E = expm(M)                               # [[exp(A Tc), integrale], [0, I]]
    PHI[j] = E[:2, :2]; GAM[j] = E[:2, 3]


def modele_pas(v, i, s, vin):
    """Un pas Tc du modele physique nominal (interpolation lineaire de la
    table en s). Renvoie v', iL', la matrice d(v', iL')/d(v, iL), les
    derivees de v' et iL' par rapport a s, puis par rapport a Vin."""
    x = s * NS
    j = min(int(x), NS - 1); f = x - j       # segment de la table et position dans le segment
    P = PHI[j] + f * (PHI[j + 1] - PHI[j]); dP = (PHI[j + 1] - PHI[j]) * NS
    g = GAM[j] + f * (GAM[j + 1] - GAM[j]); dg = (GAM[j + 1] - GAM[j]) * NS
    src = (s * vin - (1 - s) * VF) / L_MOD; dsrc = (vin + VF) / L_MOD   # source moyenne vue par la bobine
    vn = P[0, 0] * v + P[0, 1] * i + g[0] * src
    inn = P[1, 0] * v + P[1, 1] * i + g[1] * src
    dvs = dP[0, 0] * v + dP[0, 1] * i + dg[0] * src + g[0] * dsrc
    dis = dP[1, 0] * v + dP[1, 1] * i + dg[1] * src + g[1] * dsrc
    return vn, inn, P, dvs, dis, g[0] * s / L_MOD, g[1] * s / L_MOD


def pinn_pas(v, i, s, vin):
    """Un pas Tc du PINN (tau = Tc) a charge nominale, avec la matrice
    d(v', iL')/d(v, iL) et les derivees par rapport a s, calculees en mode
    direct (trois directions : v, iL, s)."""
    z = np.array([1.0, (v - 100.0) / 60.0, (i - 12.0) / 10.0, 2 * s - 1, (vin - 195.0) / 45.0,
                  (REGLAGES["G_NOM"] - G_C) / G_E])
    T = np.zeros((3, 6)); T[0, 1] = 1 / 60.0; T[1, 2] = 1 / 10.0; T[2, 3] = 2.0   # dz / d(v, iL, s)
    h = z
    for W, b in COUCHES[:-1]:
        h = np.tanh(h @ W + b)
        T = (T @ W) * (1 - h * h)
    W, b = COUCHES[-1]
    N = h @ W + b; TN = T @ W
    vn = v + SV_P * N[0]; inn = i + SI_P * N[1]
    P = np.array([[1 + SV_P * TN[0, 0], SV_P * TN[1, 0]], [SI_P * TN[0, 1], 1 + SI_P * TN[1, 1]]])
    return vn, inn, P, SV_P * TN[2, 0], SI_P * TN[2, 1], 0.0, 0.0


# Controle : le PINN relu redonne ses vecteurs de test (un pas et ses derivees)
ECART_TEST = 0.0
for q, (v0, i0, s0, vin0, G0) in enumerate(MODELE["X_TEST"]):
    vn, inn, P, dvs, dis, _, _ = pinn_pas(v0, i0, s0, vin0)
    D = np.array([[P[0, 0] - 1, P[1, 0]], [P[0, 1], P[1, 1] - 1], [dvs, dis]])
    ECART_TEST = max(ECART_TEST, abs(vn - MODELE["Y_TEST"][q, 0]), abs(inn - MODELE["Y_TEST"][q, 1]),
                     float(np.max(np.abs(D - MODELE["D_TEST"][q].reshape(3, 2)))))
print(f"Vecteurs de test du PINN : ecart maximal {ECART_TEST:.1e}.")
if ECART_TEST > 1e-12:
    raise RuntimeError("pinn_pid_modele.mat ne redonne pas ses vecteurs de test.")


class Observateur:
    """Filtre de Kalman etendu sur [Vout, iL, Vin_eff], charge nominale (inchange)."""

    def __init__(self, v0, svin):
        r = REGLAGES
        self.x = np.array([v0, 0.0, r["VIN_DEPART"]])          # depart : premiere mesure, repos, 200 V
        self.P = np.diag([r["SIG_Y"] * r["SIG_Y"], r["P0_I"], r["P0_VIN"] if svin > 0 else 0.0])
        self.Q = np.diag([r["SIG_V"] * r["SIG_V"], r["SIG_I"] * r["SIG_I"], svin * svin])
        self.R = r["SIG_Y"] * r["SIG_Y"]

    def mise_a_jour(self, y):
        """Correction par la mesure de Vout."""
        P = self.P
        S = P[0, 0] + self.R
        K = P[:, 0] / S
        self.x = self.x + K * (y - self.x[0])
        self.P = P - np.outer(K, P[0, :])

    def prediction(self, s):
        """Prediction d'un pas avec la fraction de conduction appliquee s."""
        v, i, vin = self.x
        vn, inn, Pm, _, _, dv_vin, di_vin = modele_pas(v, i, s, vin)
        F = np.array([[Pm[0, 0], Pm[0, 1], dv_vin], [Pm[1, 0], Pm[1, 1], di_vin], [0.0, 0.0, 1.0]])
        self.x = np.array([vn, inn, vin])
        self.P = F @ self.P @ F.T + self.Q


def s_continu(u, ph):
    """Fraction de conduction continue de la tranche de phase ph pour la
    commande u, et sa derivee par rapport a u (dans l'horizon)."""
    y = u * NPER - ph
    if y <= 0.0:
        return 0.0, 0.0
    if y >= NREG:
        return 1.0, 0.0
    return y / NREG, NPER / NREG


def horizon(x, x0, bloc0, r, k1, vin, pas_pred, NH, pente=None, grad=True, centre=0.0):
    """Simule NH pas a partir de l'etat x0 = (v, iL) (instant k1) avec le
    bloc PID parallele de gains K = K_ZN .* x, partant des etats du bloc
    bloc0 = (xI, xF), et le predicteur pas_pred. Renvoie le cout
      J = moyenne de 1/2 (e^2 + W_U u^2) + W_F |x - centre|^2
    (centre = 0 : sans C2 ; centre = (1, 1, 1) : correction C2) et son
    gradient par rapport a x (retropropagation dans le temps)."""
    R = REGLAGES
    pente = R["PENTE_SAT"] if pente is None else pente
    WU = R["W_U"]
    P, I, D = K_ZN * x
    xI, xF = bloc0
    v, i = x0
    pile = []                                 # valeurs gardees pour la retropropagation
    J = 0.0
    for j in range(NH):
        e = r - v                             # erreur predite
        derivee = PID_N * (D * e - xF)        # bloc PID : memes operations que BlocPID
        b = P * e + xI + derivee
        if b < D_MIN:
            u, du = D_MIN, pente
        elif b > D_MAX:
            u, du = D_MAX, pente
        else:
            u, du = b, 1.0
        entree_I = I * e
        integre = not ((b != u) and (np.sign(entree_I) == np.sign(b - u)))   # clamping
        ph = ((k1 + j) * NREG) % NPER         # phase de la porteuse
        s, ds = s_continu(u, ph)
        vn, inn, Pm, dvs, dis, _, _ = pas_pred(v, i, s, vin)
        J += 0.5 * (e * e + WU * u * u)
        pile.append((e, u, du, ds, integre, Pm, dvs, dis))
        if integre:
            xI = xI + TC * entree_I
        xF = xF + TC * derivee
        v, i = vn, inn
    dx_c = x - centre                         # ecart au centre de la regularisation
    J = J / NH + R["W_F"] * float(dx_c @ dx_c)
    if not grad:
        return J, None
    gK = np.zeros(3)                          # gradient par rapport a [P, I, D]
    lv = li = lxI = lxF = 0.0                 # adjoints de v, iL, xI, xF (etat suivant)
    for j in range(NH - 1, -1, -1):
        e, u, du, ds, integre, Pm, dvs, dis = pile[j]
        gu = WU * u / NH + (lv * dvs + li * dis) * ds
        gb = gu * du
        gder = gb + lxF * TC                  # derivee entre dans u et dans xF'
        gK[0] += gb * e
        gK[2] += gder * PID_N * e
        ge = e / NH + gb * P + gder * PID_N * D
        if integre:
            gK[1] += lxI * TC * e
            ge += lxI * TC * I
        gxI = gb + lxI
        gxF = lxF - gder * PID_N
        lv_new = Pm[0, 0] * lv + Pm[1, 0] * li - ge
        li_new = Pm[0, 1] * lv + Pm[1, 1] * li
        lv, li, lxI, lxF = lv_new, li_new, gxI, gxF
    return J, gK * K_ZN + 2 * R["W_F"] * dx_c


class BlocPID:
    """Le bloc "PID Controller" de la base commune (memes operations que
    PIDClassique du banc commun), gains P, I, D recus de l'exterieur a chaque
    pas. Conditions initiales : integrateur 0.5, filtre 0.01."""

    def __init__(self):
        self.xI = PID_CI_INTEGRATEUR
        self.xF = PID_CI_FILTRE

    def pas(self, e, K):
        P, I, D = K
        derivee = PID_N * (D * e - self.xF)
        u = P * e + self.xI + derivee
        u_sat = min(max(u, D_MIN), D_MAX)
        entree_I = I * e
        if not ((u != u_sat) and (np.sign(entree_I) == np.sign(u - u_sat))):
            self.xI += TC * entree_I
        self.xF += TC * derivee
        return u_sat


class AdaptateurPINN:
    """Observateur et optimisation des gains ; copie Python du bloc MATLAB
    pinn_pid_adaptatif.m (memes entrees, meme sortie, meme ordre des
    operations : Tester_PINN_PID_Rejeu.m compare les deux pas par pas).
    Entrees a chaque pas : erreur e, mesure de Vout, commande u sortie du
    bloc PID. Sortie : gains K = K_ZN .* x du bloc PID, qui ne dependent que
    de l'etat (changes en fin de fenetre, utilises a partir du pas suivant).
    Les etats du bloc PID (integrateur xI, filtre xF), dont l'horizon part,
    sont refaits ici a partir de e et des gains du pas, par les operations
    de BlocPID.pas (le bloc ne les donne pas a l'exterieur).
    corrections : sous-ensemble de ("C1", "C2", "C3") ; None = la
    combinaison retenue (choix_correction_pinn.json). Les autres
    arguments servent aux variantes ; leurs valeurs par defaut sont celles
    de la methode."""

    def __init__(self, predicteur="pinn", adapter=True, corrections=None, svin=None):
        r = REGLAGES
        if predicteur == "pinn":
            self.pas_pred = pinn_pas          # le PINN appris
        elif predicteur == "modele":
            self.pas_pred = modele_pas        # plafond : le modele physique nominal
        else:
            raise ValueError(predicteur)
        self.adapter = adapter
        if corrections is None:
            corrections = combinaison_retenue()
        corrections = tuple(sorted(corrections))
        if not set(corrections) <= {"C1", "C2", "C3"}:
            raise ValueError(corrections)
        self.corrections = corrections
        self.seuil = r["SEUIL_C1"] if "C1" in corrections else None   # C1 : zone morte sur l'erreur efficace
        self.centre = np.ones(3) if "C2" in corrections else np.zeros(3)  # C2 : regularisation vers ZN
        self.ensemble = "C3" in corrections                             # C3 : ensemble admissible de l'ELM-PID
        self.svin = r["SIG_VIN"] if svin is None else svin
        self.NH, self.n_adam = int(r["NH"]), int(r["N_ADAM"])
        self.x = np.ones(3)                   # multiplicateurs de Ziegler-Nichols
        self.K = K_ZN * self.x
        self.xI, self.xF = PID_CI_INTEGRATEUR, PID_CI_FILTRE   # copie des etats du bloc PID
        self.u_bloc = float("nan")            # commande recalculee par la copie du bloc
        self.obs = None
        self.k = 0                            # numero de la periode Tc
        self.cpt, self.se2 = 0, 0.0           # fenetre courante
        self.journal = []                     # une ligne par fenetre

    def mise_a_jour(self, e, mesure, u):
        k = self.k
        # 0. copie du bloc PID avec les gains de ce pas (memes operations que BlocPID.pas)
        P_, I_, D_ = self.K
        derivee = PID_N * (D_ * e - self.xF)
        b = P_ * e + self.xI + derivee
        self.u_bloc = min(max(b, D_MIN), D_MAX)
        entree_I = I_ * e
        if not ((b != self.u_bloc) and (np.sign(entree_I) == np.sign(b - self.u_bloc))):
            self.xI += TC * entree_I
        self.xF += TC * derivee
        if self.obs is None:
            self.obs = Observateur(mesure, self.svin)
        self.obs.mise_a_jour(mesure)          # 1. correction de l'observateur
        r = e + mesure                        # consigne courante
        ph = (k * NREG) % NPER                # 2. prediction de l'observateur avec la fraction appliquee
        s_app = min(max(int(np.ceil(u * NPER)) - ph, 0), NREG) / NREG
        self.obs.prediction(s_app)
        self.se2 += e * e                     # 3. fenetre
        self.cpt += 1
        self.k += 1
        if self.cpt >= NF:
            eff = np.sqrt(self.se2 / NF)
            ancien = self.x.copy()
            adapte, Jfin = False, float("nan")
            if self.adapter and (self.seuil is None or eff > self.seuil):
                Jfin = self.optimiser(r)
                adapte = True
            self.journal.append({"k": k, "eff": float(eff), "adapte": adapte, "J": Jfin,
                                 "admis": bool(g_interp(self.x) <= 0.0) if self.ensemble else None,
                                 "change": bool(np.any(self.x != ancien)), "x": self.x.copy(), "K": self.K.copy(),
                                 "vin_eff": float(self.obs.x[2]), "iL": float(self.obs.x[1])})
            self.cpt, self.se2 = 0, 0.0

    def optimiser(self, r):
        """5 iterations d'Adam sur x depuis x courant (depart a chaud),
        horizon depuis l'etat estime et l'etat du bloc, projection apres
        chaque iteration : sur la boite ou, avec C3, sur
        l'ensemble admissible de l'ELM-PID par sa projection."""
        R = REGLAGES
        x0 = (float(self.obs.x[0]), float(self.obs.x[1])); vin = float(self.obs.x[2])
        th = self.x.copy()
        m = np.zeros(3); vv = np.zeros(3)     # moments d'Adam, remis a zero (Ito et Wasa, p. 5)
        b1t, b2t = 1.0, 1.0
        J = float("nan")
        for it in range(1, self.n_adam + 1):
            b1t *= R["BETA1"]; b2t *= R["BETA2"]
            J, g = horizon(th, x0, (self.xI, self.xF), r, self.k, vin, self.pas_pred, self.NH, centre=self.centre)
            m = R["BETA1"] * m + (1 - R["BETA1"]) * g
            vv = R["BETA2"] * vv + (1 - R["BETA2"]) * g * g
            pas_adam = -R["ALPHA"] * (m / (1 - b1t)) / (np.sqrt(vv / (1 - b2t)) + R["EPS_ADAM"])
            if self.ensemble:                 # C3 : projection de l'ELM-PID (th admissible au depart)
                th = projeter(th, pas_adam, R_PROJ)
            else:                             # projection sur la boite (Algorithme 1, ligne 8)
                th = np.minimum(np.maximum(th + pas_adam, X_MIN), X_MAX)
        self.x = th
        self.K = K_ZN * th
        return float(J)


class PINNPID:
    """Montage Simulink : adaptateur -> gains -> bloc PID de la base commune
    -> u, et u renvoye a l'adaptateur (avec e et la mesure)."""
    utilise_mesure = True

    def __init__(self, **options):
        self.adapt = AdaptateurPINN(**options)
        self.bloc = BlocPID()
        self.journal = self.adapt.journal
        self.K_pas = []

    def pas(self, e, mesure):
        K = self.adapt.K                       # sortie de l'adaptateur (son etat seulement)
        self.K_pas.append(K)
        u = self.bloc.pas(e, K)
        self.adapt.mise_a_jour(e, mesure, u)
        if self.adapt.xI != self.bloc.xI or self.adapt.xF != self.bloc.xF or self.adapt.u_bloc != u:
            raise RuntimeError("La copie du bloc PID dans l'adaptateur differe du bloc PID.")
        return u


class PIDParallele(PIDClassique):
    """PID du bloc Simulink (forme Parallel, clamping) avec d'autres gains P, I, D."""

    def __init__(self, P, I, D):
        super().__init__()
        self.P, self.I, self.D = P, I, D

    def pas(self, e):
        derivee = PID_N * (self.D * e - self.xF)
        u = self.P * e + self.xI + derivee
        u_sat = min(max(u, D_MIN), D_MAX)
        entree_I = self.I * e
        if not ((u != u_sat) and (np.sign(entree_I) == np.sign(u - u_sat))):
            self.xI += TC * entree_I
        self.xF += TC * derivee
        return u_sat


# Controle 1 du gradient : hors saturation, pente de la methode
def ecart_gradient(x, x0, bloc0, r, k1, vin, pente, h_rel=1e-6, centre=0.0):
    _, g_a = horizon(x, x0, bloc0, r, k1, vin, pinn_pas, NF, pente=pente, centre=centre)
    ecart = 0.0
    for m_ in range(3):
        h_ = x[m_] * h_rel
        xp, xm = x.copy(), x.copy(); xp[m_] += h_; xm[m_] -= h_
        fd = (horizon(xp, x0, bloc0, r, k1, vin, pinn_pas, NF, pente=pente, grad=False, centre=centre)[0]
              - horizon(xm, x0, bloc0, r, k1, vin, pinn_pas, NF, pente=pente, grad=False, centre=centre)[0]) / (2 * h_)
        ecart = max(ecart, abs(fd - g_a[m_]) / abs(fd))
    return ecart


X_ESSAI = np.array([1.3, 0.7, 1.6])           # un point interieur de la boite
ECART_GRAD = ecart_gradient(X_ESSAI, (99.5, 20.0), (0.5, 7.3227e-6 * 0.5), 100.0, 7, 200.0, None)
# Controle C2 : meme gradient avec la regularisation centree sur Ziegler-Nichols, au point interieur et en un
# point proche du centre (ou le terme W_F |x - 1|^2 pese relativement plus)
ECART_GRAD_C2 = max(ecart_gradient(x_, (99.5, 20.0), (0.5, 7.3227e-6 * 0.5), 100.0, 7, 200.0, None,
                                   centre=np.ones(3)) for x_ in (X_ESSAI, np.array([1.02, 0.97, 1.03])))


def pas_en_butee(x, x0, bloc0, r, k1, vin):
    """Nombre de pas de l'horizon ou la sortie du bloc est en butee (et ou le clamping coupe l'integrateur)."""
    P_, I_, D_ = K_ZN * x
    xI, xF = bloc0; v, i = x0; n_sat = n_cl = 0
    for j in range(NF):
        e = r - v; der = PID_N * (D_ * e - xF); b = P_ * e + xI + der; u = min(max(b, D_MIN), D_MAX)
        integre = not ((b != u) and (np.sign(I_ * e) == np.sign(b - u)))
        n_sat += int(b != u); n_cl += int(not integre)
        s, _ = s_continu(u, ((k1 + j) * NREG) % NPER)
        v, i, *_ = pinn_pas(v, i, s, vin)
        if integre:
            xI += TC * I_ * e
        xF += TC * der
    return n_sat, n_cl


# Controle 2 : horizons ou la sortie passe une partie du temps en butee (haute puis basse), clamping actif,
# pente 0 (derivee exacte de la loi saturee)
CAS_SAT = [((99.0, 22.0), (0.55, 7.3227e-6)), ((102.0, 22.0), (0.5, -2 * 7.3227e-6))]
ECART_GRAD_SAT = max(ecart_gradient(X_ESSAI, x0_, b0_, 100.0, 7, 200.0, 0.0, centre=c_)
                     for x0_, b0_ in CAS_SAT for c_ in (0.0, np.ones(3)))
BUTEES = [pas_en_butee(X_ESSAI, x0_, b0_, 100.0, 7, 200.0) for x0_, b0_ in CAS_SAT]
print(f"Gradient de l'horizon contre differences finies : ecart relatif maximal {ECART_GRAD:.1e} (hors saturation), "
      f"{ECART_GRAD_SAT:.1e} (pente 0 ; pas en butee / pas sous clamping sur 110 : "
      + ", ".join(f"{a} / {b}" for a, b in BUTEES) + f") ; avec C2 (centre (1, 1, 1)) {ECART_GRAD_C2:.1e}.")
if ECART_GRAD > 1e-4 or ECART_GRAD_SAT > 1e-4 or ECART_GRAD_C2 > 1e-4:
    raise RuntimeError("Le gradient de l'horizon ne correspond pas aux differences finies.")


# Evaluation d'un regulateur sur des essais, un processus par essai
_FABRIQUE, _ESSAIS_EVAL = None, None


def _un_essai(code):
    reg = _FABRIQUE()
    sc = _ESSAIS_EVAL[code]
    sim = simuler(reg, sc)
    o = grandeurs(sim, sc)
    e = sim["consigne"] - sim["v"]
    leger = types.SimpleNamespace(journal=getattr(reg, "journal", None),
                                  K_pas=np.array(reg.K_pas) if hasattr(reg, "K_pas") else None)
    return code, {"o": o, "IAE_dem": float(np.sum(np.abs(e[:K30])) * TC),
                  "revenus": bool(all(ev["revenu"] for ev in o["evenements"])), "reg": leger, "sim": sim}


def evaluer(fabrique, codes=None, essais=None, processus=4):
    """Simule un nouveau regulateur (fabrique()) sur chaque essai."""
    global _FABRIQUE, _ESSAIS_EVAL
    _FABRIQUE, _ESSAIS_EVAL = fabrique, (essais or ESSAIS)
    codes = codes or list(_ESSAIS_EVAL)
    if processus > 1 and len(codes) > 1:
        with mp.get_context("fork").Pool(min(processus, len(codes))) as pool:
            res = dict(pool.map(_un_essai, codes, chunksize=1))
    else:
        res = dict(_un_essai(c) for c in codes)
    return {c: res[c] for c in codes}


# %% ETAPE 2 : references

def termes(res):
    """Les 13 termes du cout (10 IAE apres 30 ms, 3 IAE de demarrage)."""
    return np.array([res[c]["o"]["IAE"] for c in CODES_IAE] + [res[c]["IAE_dem"] for c in CODES_DEM])


def non_revenus(res):
    return [c for c in res if not res[c]["revenus"]]


def multiplicateurs(res, code):
    """Multiplicateurs de ZN a la fin de chaque fenetre (une ligne par fenetre)."""
    return np.array([x["x"] for x in res[code]["reg"].journal])


if __name__ == "__main__":
    titre("ETAPE 2 : les references")
    REF = {"Ziegler-Nichols": evaluer(lambda: PIDClassique()),
           "meilleur PID fige": evaluer(lambda: PIDParallele(GAINS_FIGE["P"], GAINS_FIGE["I"], GAINS_FIGE["D"]))}
    T_ZN = termes(REF["Ziegler-Nichols"])
    J_FIGE = float(np.mean(termes(REF["meilleur PID fige"]) / T_ZN))
    print(f"  Meilleur PID fige sur grille : J = {J_FIGE:.3f} (0.738 attendu).")

    # %% ETAPE 3 : PINN-PID sur les onze essais, controle d'identite

    COMB = combinaison_retenue()
    NOM_C = "PINN-PID"
    titre(f"ETAPE 3 : PINN-PID, combinaison retenue sur E1 a E4 : {'+'.join(COMB) or 'aucune'} "
          f"(sans adaptation : le PID de Ziegler-Nichols lui-meme)")
    fige = evaluer(lambda: PINNPID(adapter=False), ["S1", "S8b"])
    ECART_ZN = max(float(np.max(np.abs(fige[c]["sim"]["d"] - REF["Ziegler-Nichols"][c]["sim"]["d"]))) for c in fige)
    print(f"  Controle : PINN-PID sans adaptation contre Ziegler-Nichols (S1, S8b), ecart maximal sur d : {ECART_ZN:.1e}.")
    if ECART_ZN != 0.0:
        raise RuntimeError("Sans adaptation, le PINN-PID ne redonne pas le PID de Ziegler-Nichols a l'identique.")

    def bilan(nom, res):
        r_ = termes(res) / T_ZN
        J = float(np.mean(r_))
        non = non_revenus(res)
        print(f"\n  {nom} : J = {J:.3f} (apres 30 ms {np.mean(r_[:10]):.3f}, demarrage {np.mean(r_[10:]):.3f}) ; "
              f"{'tous les evenements reviennent' if not non else 'pas revenus : ' + ', '.join(non)}")
        for code, r in res.items():
            o = r["o"]; X = multiplicateurs(res, code)
            n_opt = sum(1 for x in r["reg"].journal[59:] if x["adapte"])
            print(f"    {code:4s} optimise apres 30 ms {n_opt:3d}/{len(r['reg'].journal) - 59} ; depassement {o['depassement_pct']:5.1f} %, IAE {o['IAE'] * 1e3:7.2f} mV.s, dem. "
                  f"{r['IAE_dem'] * 1e3:7.2f} mV.s, efficace {o['e_eff_V'] * 1e3:6.1f} mV ; x final "
                  f"({X[-1, 0]:.3f}, {X[-1, 1]:.3f}, {X[-1, 2]:.3f}), x apres 30 ms de ({X[59:, 0].min():.2f}, "
                  f"{X[59:, 1].min():.2f}, {X[59:, 2].min():.2f}) a ({X[59:, 0].max():.2f}, {X[59:, 1].max():.2f}, "
                  f"{X[59:, 2].max():.2f})", flush=True)
        return J

    t0 = time.time()
    PINN = evaluer(lambda: PINNPID())
    J_PINN = bilan(NOM_C, PINN)
    print(f"  ({time.time() - t0:.0f} s)")

    # %% ETAPE 4 : ablations (rapportees, jamais choisies) et sensibilite

    titre("ETAPE 4 : ablations (rapportees, jamais choisies)")
    X_FIN_S1 = multiplicateurs(PINN, "S1")[-1]
    ABLATIONS = {"Plafond (modele physique nominal a la place du PINN)": lambda: PINNPID(predicteur="modele"),
                 "Aucune correction (sans zone morte, 0.893 attendu)": lambda: PINNPID(corrections=()),
                 "C1+C2 (a 1 % sur E1-E4, non retenue)": lambda: PINNPID(corrections=("C1", "C2")),
                 "C1+C3 (a 1 % sur E1-E4, non retenue)": lambda: PINNPID(corrections=("C1", "C3")),
                 "PID fige aux gains finaux du PINN-PID sur S1": lambda: PIDParallele(*(K_ZN * X_FIN_S1))}
    RES_ABL, J_ABL = {}, {}
    for nom, fab in ABLATIONS.items():
        t0 = time.time()
        RES_ABL[nom] = evaluer(fab)
        r_ = termes(RES_ABL[nom]) / T_ZN
        J_ABL[nom] = float(np.mean(r_))
        print(f"  {nom:52s} J = {J_ABL[nom]:.4f} ({np.mean(r_[:10]):.3f} / {np.mean(r_[10:]):.3f}) ; non revenus : "
              f"{','.join(non_revenus(RES_ABL[nom])) or '-'} ({time.time() - t0:.0f} s)", flush=True)
    print(f"  (gains finaux S1 x({X_FIN_S1[0]:.3f}, {X_FIN_S1[1]:.3f}, {X_FIN_S1[2]:.3f}))")

    titre("ETAPE 4b : sensibilite (L ou C a +-0.1 %)")
    L_NOM, C_NOM = L_BOB, C_CONV
    SENSIBILITE = []
    for nom_p, fl, fc in (("L - 0.1 %", 0.999, 1.0), ("L + 0.1 %", 1.001, 1.0), ("C - 0.1 %", 1.0, 0.999),
                          ("C + 0.1 %", 1.0, 1.001)):
        L_BOB, C_CONV = L_NOM * fl, C_NOM * fc   # seul le circuit change ; le modele du regulateur reste nominal
        t_zn = termes(evaluer(lambda: PIDClassique()))
        r_p = evaluer(lambda: PINNPID())
        J_p = float(np.mean(termes(r_p) / t_zn))
        SENSIBILITE.append({"circuit": nom_p, "J": J_p, "non_revenus": non_revenus(r_p)})
        print(f"  {nom_p:10s} J = {J_p:.3f} ; non revenus : {','.join(non_revenus(r_p)) or '-'}", flush=True)
    L_BOB, C_CONV = L_NOM, C_NOM

    titre("ETAPE 4c : marges de phase du modele moyen aux gains finaux")
    W_F = np.logspace(np.log10(2 * np.pi * 20.0), np.log10(np.pi / TC * 0.999), 8000)
    Z_F = np.exp(1j * W_F * TC)

    def marge(K, R, V):
        """Marge de phase (degres) du bloc parallele de gains K sur le modele
        moyen sans pertes au coin (R, V) (memes formules que boite_gains_pinn.py)."""
        num, den, _ = cont2discrete(([V], [L_BOB * C_CONV, L_BOB / R, 1.0]), TC, method="zoh")
        G = np.polyval(np.squeeze(num), Z_F) / np.polyval(den, Z_F)
        Cz = K[0] + K[1] * TC / (Z_F - 1.0) + K[2] * PID_N * (Z_F - 1.0) / (Z_F - 1.0 + TC * PID_N)
        boucle = Cz * G
        lm = np.log(np.abs(boucle)); ph = np.unwrap(np.angle(boucle))
        j = np.where(np.diff(np.sign(lm)) != 0)[0]
        if len(j) == 0:
            return 180.0
        f = lm[j] / (lm[j] - lm[j + 1])
        phase = ph[j] + f * (ph[j + 1] - ph[j])
        return float(((np.degrees(phase) + 360.0) % 360.0 - 180.0).min())

    MARGES = {}
    for nom, x_ in [("depart (Ziegler-Nichols)", np.ones(3))] + [(f"final {c}", multiplicateurs(PINN, c)[-1])
                                                                for c in ("S1", "S3", "S8a", "S8b", "S9")]:
        nominal = min(marge(K_ZN * x_, R, V) for R in (4.0, 5.0, 6.0) for V in (160.0, 200.0, 240.0))
        leger = min(marge(K_ZN * x_, R, V) for R in (25.0, 98.0) for V in (160.0, 200.0, 240.0))
        MARGES[nom] = {"x": x_.tolist(), "marge_nominale_deg": nominal, "marge_legere_deg": leger}
        print(f"  {nom:26s} x({x_[0]:.3f}, {x_[1]:.3f}, {x_[2]:.3f}) : marge minimale {nominal:5.1f} deg (4 a 6 ohms), "
              f"{leger:5.1f} deg (25 et 98 ohms)")

    # %% ETAPE 5 : comparaison

    titre("ETAPE 5 : PINN-PID face aux references")
    TOUS = {"Ziegler-Nichols": REF["Ziegler-Nichols"], "meilleur PID fige": REF["meilleur PID fige"], NOM_C: PINN}
    J_TOUS = {n: float(np.mean(termes(r) / T_ZN)) for n, r in TOUS.items()}
    DECOMP = {n: (float(np.mean((termes(r) / T_ZN)[:10])), float(np.mean((termes(r) / T_ZN)[10:]))) for n, r in TOUS.items()}
    print("  " + " " * 28 + "".join(f"{n:>20s}" for n in TOUS))
    print("  " + f"{'cout J':28s}" + "".join(f"{J_TOUS[n]:20.3f}" for n in TOUS))
    print("  " + f"{'apres 30 ms / demarrage':28s}" + "".join(f"{DECOMP[n][0]:13.3f} /{DECOMP[n][1]:5.3f}" for n in TOUS))
    print("  " + f"{'evenements non revenus':28s}" + "".join(f"{(','.join(non_revenus(TOUS[n])) or '-'):>20s}" for n in TOUS))
    for code in ESSAIS:
        print("  " + f"{code + ' IAE mV.s':28s}" + "".join(f"{TOUS[n][code]['o']['IAE'] * 1e3:20.2f}" for n in TOUS))
    print("  " + f"{'S1 demarrage mV.s':28s}" + "".join(f"{TOUS[n]['S1']['IAE_dem'] * 1e3:20.2f}" for n in TOUS))

    # %% ETAPE 5b : essais complementaires, hors du cout J
    # S10 est simule ici avec le meme regulateur, pour figurer comme les onze essais dans
    # predictions_banc_pinn_pid.json (lu par Simuler_PINN_PID.m) et dans
    # banc_pinn_pid_resultats.json. Il n'entre ni dans J, ni dans les ablations ou la
    # sensibilite, calcules plus haut sur les onze essais seulement.

    ESSAIS_HORS_J = {sc["code"]: sc for sc in SCENARIOS_COMPLEMENTAIRES}
    HORS_J = {}
    if ESSAIS_HORS_J:
        titre("ETAPE 5b : essais complementaires, hors du cout J (" + ", ".join(ESSAIS_HORS_J) + ")")
        HORS_J = {"Ziegler-Nichols": evaluer(lambda: PIDClassique(), essais=ESSAIS_HORS_J),
                  NOM_C: evaluer(lambda: PINNPID(), essais=ESSAIS_HORS_J)}
        for code in ESSAIS_HORS_J:
            for n, res in HORS_J.items():
                print(f"  {code} {n:16s} {resume(res[code]['o'])}")
                for ev in res[code]["o"]["evenements"]:
                    etat = ("reste dans la bande" if ev["reste_dans_bande"] else
                            (f"retour en {ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU"))
                    print(f"      evenement a {ev['t_ms']:6.1f} ms : IAE {ev['IAE'] * 1e3:7.2f} mV.s, ecart max "
                          f"{ev['e_max_V']:6.2f} V, {etat}")
            X = multiplicateurs(HORS_J[NOM_C], code)
            print(f"  {code} {NOM_C} : x final ({X[-1, 0]:.3f}, {X[-1, 1]:.3f}, {X[-1, 2]:.3f})", flush=True)

    # %% ETAPE 6 : fichiers, resultats et figure

    titre("ETAPE 6 : fichiers")
    pas_1ms = int(round(1e-3 / TC))
    predictions = {"format": 2, "Te": TC, "corrections": list(COMB),   # format 2 des predictions
                   "regulateur": "PINN-PID (bloc PID de la base commune a gains externes, gains optimises "
                                 "par Adam, iteration en temps reel ; corrections " + ("+".join(COMB) or "aucune") +
                                 ", choisies sur E1 a E4)",
                   "reglages": REGLAGES, "boite": {"x_min": X_MIN.tolist(), "x_max": X_MAX.tolist()},
                   "K_depart": K_ZN.tolist(),
                   "controle": {"ecart_vecteurs_test": ECART_TEST, "ecart_sans_adaptation_ZN": ECART_ZN,
                                "ecart_gradient": ECART_GRAD, "ecart_gradient_saturation": ECART_GRAD_SAT,
                                "ecart_gradient_C2": ECART_GRAD_C2},
                   "essais": {}}
    for code, r in list(PINN.items()) + list(HORS_J.get(NOM_C, {}).items()):   # onze essais, puis S10
        Kp = r["reg"].K_pas
        predictions["essais"][code] = {"grandeurs": r["o"],
                                       "v_toutes_les_ms": r["sim"]["v"][::pas_1ms].tolist(),
                                       "iL_toutes_les_ms": r["sim"]["iL"][::pas_1ms].tolist(),
                                       "d_moyen_par_ms": [float(np.mean(r["sim"]["d"][i:i + pas_1ms]))
                                                          for i in range(0, len(r["sim"]["d"]) - pas_1ms + 1, pas_1ms)],
                                       "K_toutes_les_ms": Kp[::pas_1ms].tolist(),
                                       "fenetres_adaptees": sum(1 for x in r["reg"].journal if x["adapte"]),
                                       "fenetres_gains_changes": sum(1 for x in r["reg"].journal if x["change"]),
                                       "K_final": r["reg"].journal[-1]["K"].tolist()}
    with open(os.path.join(DOSSIER_PINN, "predictions_banc_pinn_pid.json"), "w", encoding="utf-8") as f:
        json.dump(predictions, f, indent=1, allow_nan=False)
    print("  predictions_banc_pinn_pid.json ecrit.")

    # Reglages du bloc MATLAB (pinn_pid_adaptatif.m) : ceux de l'adaptateur retenu, tels quels
    _A = AdaptateurPINN()
    if _A.ensemble:
        raise RuntimeError("Combinaison avec C3 : la projection de l'ELM-PID n'est pas dans pinn_pid_adaptatif.m.")
    savemat(os.path.join(DOSSIER_PINN, "pinn_pid_reglages.mat"),
            {"FORMAT": 3.0, "CORRECTIONS": "+".join(_A.corrections) or "aucune",
             "K_DEPART": K_ZN.reshape(1, 3), "X_MIN": X_MIN.reshape(1, 3), "X_MAX": X_MAX.reshape(1, 3),
             "CENTRE": _A.centre.reshape(1, 3), "SEUIL": float("nan") if _A.seuil is None else _A.seuil,
             "ENSEMBLE": 0.0, "PREDICTEUR_PINN": 1.0,
             "NH": float(_A.NH), "N_ADAM": float(_A.n_adam), "ALPHA": REGLAGES["ALPHA"],
             "BETA1": REGLAGES["BETA1"], "BETA2": REGLAGES["BETA2"], "EPS_ADAM": REGLAGES["EPS_ADAM"],
             "W_U": REGLAGES["W_U"], "W_F": REGLAGES["W_F"], "PENTE_SAT": REGLAGES["PENTE_SAT"],
             "G_NOM": REGLAGES["G_NOM"], "VIN_DEPART": REGLAGES["VIN_DEPART"], "SIG_Y": REGLAGES["SIG_Y"],
             "SIG_V": REGLAGES["SIG_V"], "SIG_I": REGLAGES["SIG_I"], "SIG_VIN": _A.svin,
             "P0_I": REGLAGES["P0_I"], "P0_VIN": REGLAGES["P0_VIN"],
             "NF": float(NF), "TC": TC, "N_FILTRE": PID_N, "D_MIN": D_MIN, "D_MAX": D_MAX,
             "CI_INTEGRATEUR": PID_CI_INTEGRATEUR, "CI_FILTRE": PID_CI_FILTRE,
             "NPER": float(NPER), "NREG": float(NREG), "L_MOD": L_MOD, "VF": VF,
             "PHI": PHI.reshape(NS + 1, 4), "GAM": GAM})
    print("  pinn_pid_reglages.mat ecrit (corrections " + ("+".join(_A.corrections) or "aucune") + ").")

    # Reference du rejeu pas a pas (Tester_PINN_PID_Rejeu.m) : entrees de l'adaptateur, gains, decisions
    rejeu = {}
    for code in CODES_REJEU:
        r = PINN[code]
        jr = r["reg"].journal
        rejeu[code] = {"e": (r["sim"]["consigne"] - r["sim"]["mesure"]).reshape(-1, 1),
                       "mesure": r["sim"]["mesure"].reshape(-1, 1), "u": r["sim"]["d"].reshape(-1, 1),
                       "K": np.array(r["reg"].K_pas),
                       "adapte": np.array([x["adapte"] for x in jr], float).reshape(-1, 1),
                       "vin_eff": np.array([x["vin_eff"] for x in jr]).reshape(-1, 1),
                       "iL_obs": np.array([x["iL"] for x in jr]).reshape(-1, 1)}
    savemat(os.path.join(DOSSIER_PINN, "reference_rejeu_pinn_pid.mat"), rejeu, do_compression=True)
    print(f"  reference_rejeu_pinn_pid.mat ecrit ({', '.join(CODES_REJEU)}).")

    def resume_essais(res):
        return {c: {"grandeurs": res[c]["o"], "IAE_dem": res[c]["IAE_dem"], "revenus": res[c]["revenus"]} for c in res}

    sortie = {"corrections": list(COMB), "reglages": REGLAGES, "boite": predictions["boite"], "J": J_TOUS, "decomposition": DECOMP,
              "J_ablations": J_ABL, "sensibilite": SENSIBILITE, "marges": MARGES,
              "references": {n: resume_essais(TOUS[n]) for n in TOUS},
              "ablations": {n: resume_essais(RES_ABL[n]) for n in RES_ABL},
              "multiplicateurs_par_fenetre": {c: multiplicateurs(PINN, c).tolist() for c in PINN},
              "essais_hors_J": {"codes": list(ESSAIS_HORS_J),
                                "note": "essais complementaires (S10), hors du cout J",
                                "references": {n: resume_essais(HORS_J[n]) for n in HORS_J},
                                "multiplicateurs_par_fenetre": {c: multiplicateurs(HORS_J[NOM_C], c).tolist()
                                                                for c in ESSAIS_HORS_J}},
              "duree_s": round(time.time() - T_DEBUT, 1)}
    with open(os.path.join(DOSSIER_PINN, "banc_pinn_pid_resultats.json"), "w", encoding="utf-8") as f:
        json.dump(sortie, f, indent=1, allow_nan=False, default=lambda x: None)
    print("  banc_pinn_pid_resultats.json ecrit.")

    fig, axes = plt.subplots(3, 4, figsize=(20, 9.5), sharex="col")
    for col, (code, xlim) in enumerate((("S1", (0, 30)), ("S2", (40, 90)), ("S3", (40, 90)), ("S8b", None))):
        for nom, couleur in (("Ziegler-Nichols", "0.6"), ("meilleur PID fige", "tab:green"), (NOM_C, "tab:blue")):
            sim = TOUS[nom][code]["sim"]
            axes[0, col].plot(sim["t"] * 1e3, sim["v"], lw=0.6, color=couleur, label=nom)
            axes[2, col].plot(sim["t"] * 1e3, sim["iL"], lw=0.4, color=couleur)
        axes[0, col].plot(sim["t"] * 1e3, sim["consigne"], "k--", lw=0.6)
        Kx = PINN[code]["reg"].K_pas / K_ZN
        t_ms = PINN[code]["sim"]["t"] * 1e3
        for i, nom in enumerate(("P", "I", "D")):
            axes[1, col].plot(t_ms, Kx[:, i], lw=1.0, label=nom)
        axes[0, col].set_ylim((85, 115) if code in ("S2", "S3") else (30, 200))
        axes[0, col].set_title(f"{code} : tension de sortie", fontsize=9, loc="left")
        axes[1, col].set_title(f"{code} : gains du PINN-PID (multiples de ZN)", fontsize=9, loc="left")
        axes[2, col].set_title(f"{code} : courant de la bobine", fontsize=9, loc="left")
        axes[2, col].set_xlabel("temps (ms)")
        if xlim:
            axes[2, col].set_xlim(*xlim)
        for ax in axes[:, col]:
            ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7)
    axes[1, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(DOSSIER_PINN, "banc_pinn_pid.png"), dpi=110)
    plt.close(fig)
    print(f"  banc_pinn_pid.png ecrit. Duree totale : {time.time() - T_DEBUT:.0f} s.")
