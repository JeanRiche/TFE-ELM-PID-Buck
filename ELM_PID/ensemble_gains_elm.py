# =============================================================================
# ensemble_gains_elm.py
#
# OBJECTIF
#   Construire, une fois pour toutes et hors ligne, l'ensemble des gains
#   admissibles de l'ELM-PID, dans lequel la loi d'adaptation projette les
#   gains apres chaque pas.
#
# LE CRITERE
#   Modele moyen du Buck, discretise avec un bloqueur d'ordre zero a Tc,
#   G(s) = Vin / (L C s^2 + (L/R) s + 1), boucle C(z) G(z) avec le bloc PID
#   de la base commune (forme Parallel, Forward Euler a Tc pour
#   l'integrateur et le filtre de derivee) :
#     C(z) = P + I Tc/(z-1) + D N (z-1)/(z-1+Tc N).
#   Un reglage est admissible si, aux neuf coins nominaux (R = 4, 5, 6 ohms ;
#   Vin = 160, 200, 240 V) :
#     1. sa plus petite marge de phase est au moins celle du point de depart
#        (gains de Ziegler-Nichols dans ce bloc ; la valeur est imprimee) ;
#     2. sa plus haute frequence de coupure est au plus fs/10 = 2.2 kHz.
#   La regle interdit a l'adaptation de reduire la marge nominale sous celle
#   du regulateur dont elle part ; elle ne porte que sur la plage nominale :
#   l'adaptation doit trouver seule les gains de la charge legere.
#
# PROJECTION SUR L'ENSEMBLE ADMISSIBLE
#   L'ensemble admissible lui-meme est tabule (et non le plus grand pave a axes
#   paralleles qui le contient). Comme le depart est sur la frontiere (sa
#   marge est le seuil), un pave l'aurait pour coin : Ki ne pourrait que
#   baisser, Kp et Kd que monter. La loi d'adaptation projette donc sur
#   l'ensemble en glissant le long de sa frontiere (projection des parametres
#   sur un ensemble admissible : Goodwin et Sin 1984 ; Ioannou et Sun 1996),
#   avec restauration vers la frontiere quand le pas tangent en sort (Rosen
#   1961 ; code dans banc_elm_pid.py et elm_pid_adaptatif.m).
#
# LA TABLE
#   Multiplicateurs (a, b, c) des gains de Ziegler-Nichols, chacun de 1/4 a
#   4 (deux octaves de part et d'autre, comme la grille d'origine), grille
#   de pas 2^(1/8) en echelle logarithmique (33 points par axe). En chaque
#   point : g = max( (M0 - marge) / M0 , (fc - fc_max) / fc_max ), M0 = marge
#   du depart. Admissible <=> g <= 0. Entre les points, g est interpole
#   (trilineaire en log2 des multiplicateurs) par la loi d'adaptation.
#
# CE QUE PRODUIT CE SCRIPT
#   ensemble_gains_elm.mat (table lue par le banc et le bloc MATLAB),
#   ensemble_gains_elm.json (resume et controles).
#
# COMMENT LANCER CE SCRIPT
#   python ensemble_gains_elm.py      (une a deux minutes)
# =============================================================================


# %% ETAPE 0 : banc commun, modele moyen

import os
import json
import time
import numpy as np
from scipy.signal import cont2discrete
from scipy.io import savemat

try:
    DOSSIER_ELM = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_ELM = os.getcwd()
with open(os.path.join(DOSSIER_ELM, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])

T0 = time.time()
K_ZN = np.array([PID_P, PID_I, PID_D])   # depart : gains de Ziegler-Nichols du bloc PID [P, I, D]
FC_MAXI = FSW / 10.0
COINS_NOMINAUX = [(R, V) for R in (4.0, 5.0, 6.0) for V in (160.0, 200.0, 240.0)]
COINS_LEGERS = [(R, V) for R in (25.0, 98.0) for V in (160.0, 200.0, 240.0)]
EXPOSANTS = np.arange(-16, 17) / 8.0          # log2 des multiplicateurs : -2 a 2 par pas de 1/8
W = np.logspace(np.log10(2 * np.pi * 20.0), np.log10(np.pi / TC * 0.999), 8000)
Z = np.exp(1j * W * TC)


def reponse_buck(R, V):
    num_d, den_d, _ = cont2discrete(([V], [L_BOB * C_CONV, L_BOB / R, 1.0]), TC, method="zoh")
    return np.polyval(np.squeeze(num_d), Z) / np.polyval(den_d, Z)


def preparer(coins):
    G = np.array([reponse_buck(R, V) for R, V in coins])
    return (K_ZN[0] * G, K_ZN[1] * TC / (Z - 1.0) * G,
            K_ZN[2] * PID_N * (Z - 1.0) / (Z - 1.0 + TC * PID_N) * G)


TERMES_NOM, TERMES_LEG = preparer(COINS_NOMINAUX), preparer(COINS_LEGERS)


def frequentiel(a, b, c, termes):
    """Plus petite marge de phase (degres) et plus haute coupure (Hz) sur les
    coins."""
    tp, ti, td = termes
    boucle = a * tp + b * ti + c * td
    lm = np.log(np.abs(boucle))
    ph = np.unwrap(np.angle(boucle), axis=1)
    lignes, j = np.where(np.diff(np.sign(lm), axis=1) != 0)
    if len(j) == 0:
        return 180.0, 0.0
    f = lm[lignes, j] / (lm[lignes, j] - lm[lignes, j + 1])
    phase = ph[lignes, j] + f * (ph[lignes, j + 1] - ph[lignes, j])
    wc = W[j] * (W[j + 1] / W[j]) ** f
    marges = (np.degrees(phase) + 360.0) % 360.0 - 180.0
    return float(marges.min()), float(wc.max() / (2.0 * np.pi))


M0, FC0 = frequentiel(1.0, 1.0, 1.0, TERMES_NOM)
print(f"Depart (Ziegler-Nichols, bloc PID parallele) : marge nominale minimale {M0:.3f} degres, coupure {FC0:.0f} Hz.")


def g_de(a, b, c):
    m, fc = frequentiel(a, b, c, TERMES_NOM)
    return max((M0 - m) / M0, (fc - FC_MAXI) / FC_MAXI), m, fc


# %% ETAPE 1 : la table

n = len(EXPOSANTS)
G = np.zeros((n, n, n))
MARGE = np.zeros((n, n, n))
FC = np.zeros((n, n, n))
for i, ea in enumerate(EXPOSANTS):
    for j, eb in enumerate(EXPOSANTS):
        for k, ec in enumerate(EXPOSANTS):
            G[i, j, k], MARGE[i, j, k], FC[i, j, k] = g_de(2.0 ** ea, 2.0 ** eb, 2.0 ** ec)
i0 = int(np.where(EXPOSANTS == 0)[0][0])
G[i0, i0, i0] = 0.0                           # depart exactement sur la frontiere (g = 0 au lieu de +-1e-16)
print(f"Table : {n ** 3} reglages, {int(np.sum(G <= 0))} admissibles ({np.mean(G <= 0):.1%}) ({time.time() - T0:.0f} s).")


def g_interp(x):
    """Interpolation trilineaire de G en log2 des multiplicateurs x = (a, b, c)
    (meme calcul dans le bloc MATLAB)."""
    u = (np.log2(np.asarray(x, float)) - EXPOSANTS[0]) * 8.0
    u = np.minimum(np.maximum(u, 0.0), n - 1.0)
    i = np.minimum(np.floor(u).astype(int), n - 2)
    f = u - i
    s = 0.0
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = (f[0] if di else 1 - f[0]) * (f[1] if dj else 1 - f[1]) * (f[2] if dk else 1 - f[2])
                s += w * G[i[0] + di, i[1] + dj, i[2] + dk]
    return float(s)


# %% ETAPE 2 : controles

# 1. Entre les points de la table : points tires au hasard, g interpole et g vrai
rng = np.random.default_rng(1)
x = 2.0 ** rng.uniform(-2, 2, (3000, 3))
gi = np.array([g_interp(p) for p in x])
gv = np.array([g_de(*p)[0] for p in x])
faux_admis = int(np.sum((gi <= 0) & (gv > 0)))
pire = float(np.max(np.where(gi <= 0, gv, -np.inf)))
print(f"Controle 1 (3000 points hors grille) : {faux_admis} admis a tort par l'interpolation ; g vrai au pire "
      f"{pire:.4f} parmi les admis (marge {M0 * (1 - pire):.2f} degres si c'est la marge qui limite).")
# 2. Directions depuis le depart
print("Controle 2 : depuis le depart, g vrai a +-10 % sur chaque gain (negatif = admissible) :")
dirs = {}
for nom, v in (("Kp", 0), ("Ki", 1), ("Kd", 2)):
    for s in (0.9, 1.1):
        p = np.ones(3)
        p[v] = s
        dirs[f"{nom} x{s}"] = g_de(*p)[0]
        print(f"    {nom} x {s} : g = {dirs[f'{nom} x{s}']:+.4f}")
p = np.array([1.0, 1.2, 1.3])
print(f"    Ki x 1.2 avec Kd x 1.3 : g = {g_de(*p)[0]:+.4f} (mouvement couple que la boite d'origine interdisait)")


# %% ETAPE 3 : ecriture

savemat(os.path.join(DOSSIER_ELM, "ensemble_gains_elm.mat"),
        {"G_TABLE": G, "LOG2_MIN": float(EXPOSANTS[0]), "PAS_LOG2": 1.0 / 8.0, "N_TABLE": float(n),
         "MULT_MIN": 2.0 ** EXPOSANTS[0], "MULT_MAX": 2.0 ** EXPOSANTS[-1], "K_DEPART": K_ZN.reshape(1, -1),
         "MARGE_DEPART": M0, "FC_MAXI": FC_MAXI, "FORMAT": 2.0})
with open(os.path.join(DOSSIER_ELM, "ensemble_gains_elm.json"), "w", encoding="utf-8") as f:
    json.dump({"format": 2, "loi": "bloc PID parallele de la base commune", "K_depart": K_ZN.tolist(), "marge_depart_deg": M0, "coupure_depart_Hz": FC0,
               "fc_maxi_Hz": FC_MAXI, "coins_nominaux": COINS_NOMINAUX, "exposants_log2": EXPOSANTS.tolist(),
               "part_admissible": float(np.mean(G <= 0)), "controle_hors_grille": {"faux_admis": faux_admis,
                                                                                   "g_vrai_pire_admis": pire},
               "directions_depuis_depart": dirs, "duree_s": round(time.time() - T0, 1)}, f, indent=1)
print(f"ensemble_gains_elm.mat et ensemble_gains_elm.json ecrits ({time.time() - T0:.0f} s).")
