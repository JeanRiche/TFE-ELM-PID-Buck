# =============================================================================
# stabilite.py
#
# OBJECTIF
#   Analyse de stabilite lineaire de la boucle a gains figes : convertisseur
#   seul, PID de Ziegler-Nichols, ELM-PID sur tout son ensemble admissible.
#   Les gains et l'ensemble admissible sont lus tels quels.
#
# MODELE
#   Modele moyen sans pertes du Buck, etat [iL ; vC], discretise par bloqueur
#   d'ordre zero a Tc ; aucun retard de calcul (le banc applique d(k) a
#   l'instant ou v(k Tc) est lu) ; bloc PID commun
#     C(z) = P + I Tc/(z-1) + D N (z-1)/(z-1+N Tc).
#   Variante de sensibilite : un retard d'un pas Tc (z^-1).
#
# PRODUIT
#   resultats_stabilite.md, resultats_stabilite.json, figures/*.pdf
#
# LANCER : python stabilite.py   (quelques minutes)
# =============================================================================

import os
import json
import time
import numpy as np
from scipy.signal import cont2discrete, ss2tf
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

T0 = time.time()
ICI = os.path.dirname(os.path.abspath(__file__))
RACINE = os.path.dirname(ICI)
DOSSIER_ELM = os.path.join(RACINE, "ELM_PID")
FIG = os.path.join(ICI, "figures")
os.makedirs(FIG, exist_ok=True)

# %% ETAPE 0 : constantes du banc commun (lues dans banc_commun.py, pas recopiees)
with open(os.path.join(DOSSIER_ELM, "banc_commun.py"), encoding="utf-8") as f:
    SRC = f.read()
ESPACE = {"__file__": os.path.join(DOSSIER_ELM, "banc_commun.py"), "__name__": "banc"}
exec(SRC[:SRC.index("# === FIN DES DEFINITIONS DU BANC COMMUN ===")], ESPACE)
L_BOB, C_CONV, TC, FSW = ESPACE["L_BOB"], ESPACE["C_CONV"], ESPACE["TC"], ESPACE["FSW"]
P0, I0, D0, NF = ESPACE["PID_P"], ESPACE["PID_I"], ESPACE["PID_D"], ESPACE["PID_N"]
K_ZN = np.array([P0, I0, D0])
ENS = loadmat(os.path.join(DOSSIER_ELM, "ensemble_gains_elm.mat"))
with open(os.path.join(DOSSIER_ELM, "ensemble_gains_elm.json"), encoding="utf-8") as f:
    ENS_JSON = json.load(f)
G_TABLE = ENS["G_TABLE"]
LOG2_MIN, PAS_LOG2 = float(ENS["LOG2_MIN"].item()), float(ENS["PAS_LOG2"].item())
NT = int(ENS["N_TABLE"].item())
M0 = float(ENS["MARGE_DEPART"].item())
FC_MAXI = float(ENS["FC_MAXI"].item())
assert np.allclose(ENS["K_DEPART"].ravel(), K_ZN), "gains de depart differents du banc"

COINS = [(R, V) for R in (4.0, 5.0, 6.0) for V in (160.0, 200.0, 240.0)]
COINS_LEGERS = [(R, V) for R in (25.0, 98.0) for V in (160.0, 200.0, 240.0)]
NOMINAL = (5.0, 200.0)
A_F = 1.0 - NF * TC                                   # pole du filtre de derivee
W = np.logspace(np.log10(2 * np.pi * 20.0), np.log10(np.pi / TC * 0.999), 8000)   # meme grille que ensemble_gains_elm.py
ZW = np.exp(1j * W * TC)
RES = {}


def virg(x, nd=2):
    """Nombre au format francais (virgule decimale, espace des milliers)."""
    s = f"{x:,.{nd}f}".replace(",", " ").replace(".", ",")
    return s.replace("-", "\u2212")


# %% ETAPE 1 : convertisseur seul

def etat(R, V):
    A = np.array([[0.0, -1.0 / L_BOB], [1.0 / C_CONV, -1.0 / (R * C_CONV)]])
    B = np.array([[V / L_BOB], [0.0]])
    C = np.array([[0.0, 1.0]])
    D = np.array([[0.0]])
    return A, B, C, D


def procede_discret(R, V):
    A, B, C, D = etat(R, V)
    Ad, Bd, Cd, Dd, _ = cont2discrete((A, B, C, D), TC, method="zoh")
    num, den = ss2tf(Ad, Bd, Cd, Dd)
    return np.real(num.ravel()), np.real(den.ravel()), Ad


seul = []
for R, V in COINS + COINS_LEGERS:
    A, B, C, D = etat(R, V)
    p = np.sort(np.linalg.eigvals(A).astype(complex))
    wn = 1.0 / np.sqrt(L_BOB * C_CONV)
    zeta = np.sqrt(L_BOB / C_CONV) / (2.0 * R)
    _, den_d, Ad = procede_discret(R, V)
    pz = np.linalg.eigvals(Ad)
    ecart_zoh = float(np.max(np.abs(np.sort_complex(pz) - np.sort_complex(np.exp(p * TC)))))
    seul.append({"R": R, "Vin": V, "poles_s": [[float(x.real), float(x.imag)] for x in p], "zeta": zeta,
                 "wn": wn, "poles_z_module": [float(abs(x)) for x in pz], "ecart_zoh": ecart_zoh})
A5, B5, C5, D5 = etat(*NOMINAL)
p5 = np.sort(np.linalg.eigvals(A5).real)
RES["convertisseur"] = {"A": A5.tolist(), "B": B5.tolist(), "C": C5.tolist(), "D": D5.tolist(),
                        "poles_nominal": p5.tolist(), "coins": seul}
print(f"Convertisseur 5 ohms : poles {p5[1]:.2f} et {p5[0]:.2f} s^-1, zeta {seul[4]['zeta']:.3f}, wn {seul[4]['wn']:.1f} rad/s.")


# %% ETAPE 2 : boucle fermee, outils vectorises

def polys_boucle(R, V, retard=False):
    """Polynome caracteristique de la boucle fermee, lineaire en x :
    chi(z) = base + xP q1 + xI q2 + xD q3 (coefficients, puissances decroissantes)."""
    numG, denG, _ = procede_discret(R, V)
    numG = np.trim_zeros(numG, "f")
    if retard:                                        # z^-1 sur la commande : G(z)/z
        denG = np.polymul(denG, [1.0, 0.0])
    denC = np.polymul([1.0, -1.0], [1.0, -A_F])
    p1 = P0 * denC
    p2 = I0 * TC * np.array([0.0, 1.0, -A_F])
    p3 = D0 * NF * np.polymul([1.0, -1.0], [1.0, -1.0])
    base = np.polymul(denC, denG)
    n = len(base)

    def cale(q):
        q = np.polymul(q, numG)
        return np.concatenate([np.zeros(n - len(q)), q])
    return base, cale(p1), cale(p2), cale(p3)


def modules_max(X, R, V, retard=False):
    """Module maximal des poles de la boucle fermee pour chaque ligne de X (m, 3)."""
    base, q1, q2, q3 = polys_boucle(R, V, retard)
    co = base[None, :] + X[:, :1] * q1 + X[:, 1:2] * q2 + X[:, 2:3] * q3
    co = co / co[:, :1]
    deg = co.shape[1] - 1
    out = np.empty(len(X))
    for d0 in range(0, len(X), 200000):
        c = co[d0:d0 + 200000]
        M = np.zeros((len(c), deg, deg))
        M[:, 0, :] = -c[:, 1:]
        M[:, np.arange(1, deg), np.arange(0, deg - 1)] = 1.0
        out[d0:d0 + 200000] = np.max(np.abs(np.linalg.eigvals(M)), axis=1)
    return out


def poles(x, R, V, retard=False):
    base, q1, q2, q3 = polys_boucle(R, V, retard)
    return np.roots(base + x[0] * q1 + x[1] * q2 + x[2] * q3)


def termes(R, V, z=ZW, retard=False):
    numG, denG, _ = procede_discret(R, V)
    G = np.polyval(numG, z) / np.polyval(denG, z)
    if retard:
        G = G / z
    return np.array([P0 * G, I0 * TC / (z - 1.0) * G, D0 * NF * (z - 1.0) / (z - 1.0 + TC * NF) * G])


def marges(X, T, Tny):
    """Pour X (m, 3) et un coin (termes T (3, nw), termes a Nyquist Tny (3,)) :
    marge de phase min (deg), coupure max (Hz), marge de gain min (dB) et
    sa frequence (Hz), marge de module min |1+L|, marge de gain a la baisse (dB). Meme calcul de la marge
    de phase que ensemble_gains_elm.py."""
    Lw = X @ T                                        # (m, nw)
    m = len(X)
    lm = np.log(np.abs(Lw))
    ph = np.angle(Lw)
    MP = np.full(m, 180.0)
    FCm = np.zeros(m)
    li, j = np.where(np.diff(np.sign(lm), axis=1) != 0)
    if len(j):
        f = lm[li, j] / (lm[li, j] - lm[li, j + 1])
        phu = np.unwrap(np.stack([ph[li, j], ph[li, j + 1]], axis=1), axis=1)
        phase = phu[:, 0] + f * (phu[:, 1] - phu[:, 0])
        wc = W[j] * (W[j + 1] / W[j]) ** f
        mg = (np.degrees(phase) + 360.0) % 360.0 - 180.0
        np.minimum.at(MP, li, mg)
        np.maximum.at(FCm, li, wc / (2 * np.pi))
    # marges de gain : franchissements de l'axe reel negatif. |L| < 1 : marge a la hausse
    # (facteur dont le gain peut monter) ; |L| > 1 : marge a la baisse (boucle conditionnellement
    # stable, facteur dont le gain peut descendre), notee positive en dB.
    GM = np.full(m, np.inf)
    WG = np.full(m, np.nan)
    GB = np.full(m, np.inf)
    im, re = Lw.imag, Lw.real
    li, j = np.where(np.diff(np.sign(im), axis=1) != 0)
    if len(j):
        f = im[li, j] / (im[li, j] - im[li, j + 1])
        r = re[li, j] + f * (re[li, j + 1] - re[li, j])
        a = np.abs(Lw[li, j]) * (np.abs(Lw[li, j + 1]) / np.abs(Lw[li, j])) ** f
        g = -20 * np.log10(a)
        wg = W[j] * (W[j + 1] / W[j]) ** f
        for k_ in np.where(r < 0)[0]:
            if g[k_] >= 0 and g[k_] < GM[li[k_]]:
                GM[li[k_]], WG[li[k_]] = g[k_], wg[k_] / (2 * np.pi)
            elif g[k_] < 0:
                GB[li[k_]] = min(GB[li[k_]], -g[k_])
    Lny = X @ Tny                                     # z = -1 : L reel
    gn = -20 * np.log10(np.abs(Lny))
    sel = (Lny.real < 0) & (gn >= 0) & (gn < GM)
    GM[sel], WG[sel] = gn[sel], FSW * 5.0
    MM = np.min(np.abs(1.0 + Lw), axis=1)
    return MP, FCm, GM, WG, MM, GB


TERMES = {c: termes(*c) for c in COINS + COINS_LEGERS}
TERMES_NY = {c: termes(*c, z=np.array([-1.0 + 0j]))[:, 0] for c in COINS + COINS_LEGERS}
TERMES_R = {c: termes(*c, retard=True) for c in COINS}
TERMES_NY_R = {c: termes(*c, z=np.array([-1.0 + 0j]), retard=True)[:, 0] for c in COINS}


# %% ETAPE 3 : Ziegler-Nichols

UN = np.ones((1, 3))
zn = []
for c in COINS + COINS_LEGERS:
    mp, fc, gm, wg, mm, gb = marges(UN, TERMES[c], TERMES_NY[c])
    pz = poles(np.ones(3), *c)
    ligne = {"R": c[0], "Vin": c[1], "module_max": float(np.max(np.abs(pz))),
             "poles": [[float(z.real), float(z.imag)] for z in pz],
             "marge_phase_deg": float(mp[0]), "coupure_Hz": float(fc[0]), "marge_gain_dB": float(gm[0]),
             "f_180_Hz": float(wg[0]), "marge_module": float(mm[0]),
             "marge_gain_baisse_dB": float(gb[0])}
    if c in COINS:
        mpr, fcr, gmr, wgr, _, _ = marges(UN, TERMES_R[c], TERMES_NY_R[c])
        ligne.update({"retard_marge_phase_deg": float(mpr[0]), "retard_marge_gain_dB": float(gmr[0]),
                      "retard_module_max": float(np.max(np.abs(poles(np.ones(3), *c, retard=True))))})
    zn.append(ligne)
RES["ziegler_nichols"] = zn


def gain_critique(x, R, V, haut=True):
    """Facteur k (k > 1 a la hausse, k < 1 a la baisse) qui amene un pole sur le
    cercle unite, cherche sur les poles (verification independante des marges de gain)."""
    ks = np.logspace(0, 4, 2000) if haut else np.logspace(0, -4, 2000)
    mods = modules_max(ks[:, None] * np.asarray(x)[None, :], R, V)
    i = np.where(mods >= 1.0)[0]
    if len(i) == 0 or i[0] == 0:
        return float("nan")
    a, b = ks[i[0] - 1], ks[i[0]]
    for _ in range(60):
        c = np.sqrt(a * b)
        if modules_max(np.array([c * np.asarray(x)]), R, V)[0] >= 1.0:
            b = c
        else:
            a = c
    return float(20 * np.log10(np.sqrt(a * b)))


for z in zn:
    z["marge_gain_poles_dB"] = gain_critique(np.ones(3), z["R"], z["Vin"])
    if z["module_max"] >= 1.0:
        z["marge_gain_baisse_poles_dB"] = float("nan")
    else:
        z["marge_gain_baisse_poles_dB"] = -gain_critique(np.ones(3), z["R"], z["Vin"], haut=False)
znn = [z for z in zn if (z["R"], z["Vin"]) == NOMINAL][0]
print(f"ZN nominal : marge {znn['marge_phase_deg']:.2f} deg a {znn['coupure_Hz']:.0f} Hz, marge de gain "
      f"{znn['marge_gain_dB']:.2f} dB a {znn['f_180_Hz']:.0f} Hz, module max {znn['module_max']:.5f}.")


# %% ETAPE 4 : ELM-PID, ensemble admissible

def g_interp_vect(X):
    """g_interp de banc_elm_pid.py, vectorise (memes operations)."""
    u = (np.log2(X) - LOG2_MIN) / PAS_LOG2
    u = np.minimum(np.maximum(u, 0.0), NT - 1.0)
    i = np.minimum(np.floor(u), NT - 2.0).astype(int)
    f = u - i
    s = np.zeros(len(X))
    for di in (0, 1):
        wi = f[:, 0] if di else 1.0 - f[:, 0]
        for dj in (0, 1):
            wj = f[:, 1] if dj else 1.0 - f[:, 1]
            for dk in (0, 1):
                wk = f[:, 2] if dk else 1.0 - f[:, 2]
                s += wi * wj * wk * G_TABLE[i[:, 0] + di, i[:, 1] + dj, i[:, 2] + dk]
    return s


def grille(pas_par_octave):
    e = np.arange(-2 * pas_par_octave, 2 * pas_par_octave + 1) / pas_par_octave
    E = np.stack(np.meshgrid(e, e, e, indexing="ij"), axis=-1).reshape(-1, 3)
    X = 2.0 ** E
    return X, g_interp_vect(X)


def pire_marges(X, termes_d=TERMES, ny_d=TERMES_NY, coins=COINS, bloc=40):
    m = len(X)
    out = {k: np.full(m, np.inf) for k in ("MP", "GM", "MM", "GB")}
    out["FC"] = np.zeros(m)
    out["coin_MP"] = np.zeros(m, int)
    for d0 in range(0, m, bloc):
        Xc = X[d0:d0 + bloc]
        for ic, c in enumerate(coins):
            mp, fc, gm, _, mm, gb = marges(Xc, termes_d[c], ny_d[c])
            s = slice(d0, d0 + len(Xc))
            pire = mp < out["MP"][s]
            out["coin_MP"][s] = np.where(pire, ic, out["coin_MP"][s])
            out["MP"][s] = np.minimum(out["MP"][s], mp)
            out["GM"][s] = np.minimum(out["GM"][s], gm)
            out["MM"][s] = np.minimum(out["MM"][s], mm)
            out["GB"][s] = np.minimum(out["GB"][s], gb)
            out["FC"][s] = np.maximum(out["FC"][s], fc)
    return out


def pire_modules(X, coins=COINS, retard=False):
    M = np.zeros(len(X))
    arg = np.zeros(len(X), int)
    for ic, c in enumerate(coins):
        mod = modules_max(X, *c, retard=retard)
        arg = np.where(mod > M, ic, arg)
        M = np.maximum(M, mod)
    return M, arg


def fmt_x(x):
    return "(" + " ; ".join(virg(v, 3) for v in x) + ")"


elm = {}
# (i) points admissibles de la table (pas 2^(1/8))
X8, g8 = grille(8)
A8 = X8[g8 <= 0]
t = time.time()
M8 = pire_marges(A8)
mod8, _ = pire_modules(A8)
elm["table"] = {"n": int(len(A8)), "pire_MP": float(M8["MP"].min()), "x_pire_MP": A8[np.argmin(M8["MP"])].tolist(),
                "coin_pire_MP": COINS[int(M8["coin_MP"][np.argmin(M8["MP"])])],
                "pire_GM": float(M8["GM"].min()), "x_pire_GM": A8[np.argmin(M8["GM"])].tolist(),
                "fc_max": float(M8["FC"].max()), "pire_MM": float(M8["MM"].min()),
                "pire_GB": float(M8["GB"].min()), "n_conditionnels": int(np.sum(np.isfinite(M8["GB"]))),
                "pire_module": float(mod8.max()), "x_pire_module": A8[np.argmax(mod8)].tolist(),
                "n_instables": int(np.sum(mod8 >= 1.0))}
print(f"Table : {len(A8)} points admissibles, pire marge {M8['MP'].min():.3f} deg, pire module {mod8.max():.6f} ({time.time()-t:.0f} s).")

# (ii) grille fine 2^(1/32) pour les poles
t = time.time()
X32, g32 = grille(32)
A32 = X32[g32 <= 0]
mod32, arg32 = pire_modules(A32)
i_pm = int(np.argmax(mod32))
elm["fine_poles"] = {"n_grille": int(len(X32)), "n_admis": int(len(A32)), "pire_module": float(mod32.max()),
                     "x_pire_module": A32[i_pm].tolist(), "coin_pire_module": COINS[int(arg32[i_pm])],
                     "n_instables": int(np.sum(mod32 >= 1.0)),
                     "quantiles_module": [float(q) for q in np.quantile(mod32, [0.5, 0.9, 0.99])]}
print(f"Grille fine 2^(1/32) : {len(A32)} admis, pire module {mod32.max():.6f}, instables {np.sum(mod32 >= 1)} ({time.time()-t:.0f} s).")

# (iii) grille 2^(1/16) pour les marges
t = time.time()
X16, g16 = grille(16)
A16 = X16[g16 <= 0]
M16 = pire_marges(A16)
mod16, _ = pire_modules(A16)
i_mp = int(np.argmin(M16["MP"]))
i_gm = int(np.argmin(M16["GM"]))
elm["fine_marges"] = {"n_admis": int(len(A16)), "pire_MP": float(M16["MP"].min()), "x_pire_MP": A16[i_mp].tolist(),
                      "coin_pire_MP": COINS[int(M16["coin_MP"][i_mp])],
                      "n_sous_M0": int(np.sum(M16["MP"] < M0 - 1e-9)),
                      "pire_GM": float(M16["GM"].min()), "x_pire_GM": A16[i_gm].tolist(),
                      "fc_max": float(M16["FC"].max()), "pire_MM": float(M16["MM"].min()),
                      "pire_GB": float(M16["GB"].min()), "x_pire_GB": A16[int(np.argmin(M16["GB"]))].tolist(),
                      "n_conditionnels": int(np.sum(np.isfinite(M16["GB"]))),
                      "x_sous_M0": A16[M16["MP"] < M0 - 1e-9].tolist(),
                      "g_sous_M0": g16[g16 <= 0][M16["MP"] < M0 - 1e-9].tolist(),
                      "n_GM_sous_6dB": int(np.sum(M16["GM"] < 6.0)),
                      "part_GM_sous_ZN": float(np.mean(M16["GM"] < min(z["marge_gain_dB"] for z in zn[:9]) - 1e-9))}
print(f"Grille 2^(1/16) : {len(A16)} admis, pire marge {M16['MP'].min():.3f} deg, pire marge de gain "
      f"{M16['GM'].min():.2f} dB ({time.time()-t:.0f} s).")

# (iv) frontiere : points admis de la grille 2^(1/16) dont un voisin ne l'est pas
n16 = 65
adm = (g16 <= 0).reshape(n16, n16, n16)
fr = np.zeros_like(adm)
for ax in range(3):
    for s in (1, -1):
        v = np.roll(adm, s, axis=ax)
        bord = np.zeros_like(adm)
        idx = [slice(None)] * 3
        idx[ax] = 0 if s == 1 else -1
        bord[tuple(idx)] = True
        fr |= adm & (~v) & (~bord)
fr_adm = fr.reshape(-1)[g16 <= 0]
elm["frontiere"] = {"n": int(fr_adm.sum()), "pire_MP_frontiere": float(M16["MP"][fr_adm].min()),
                    "pire_MP_interieur": float(M16["MP"][~fr_adm].min()),
                    "pire_module_frontiere": float(mod16[fr_adm].max()),
                    "pire_module_interieur": float(mod16[~fr_adm].max())}

# (v) gains reellement atteints par l'ELM-PID
with open(os.path.join(DOSSIER_ELM, "banc_elm_pid_resultats.json"), encoding="utf-8") as f:
    BANC = json.load(f)
vus = []
for code, liste in list(BANC["gains_elm"].items()) + list(BANC["essais_hors_J"]["gains_elm"].items()):
    for k in liste:
        vus.append(np.array(k) / K_ZN)
XV = np.unique(np.round(np.array(vus), 12), axis=0)
MV = pire_marges(XV)
modv, _ = pire_modules(XV)
gv = g_interp_vect(XV)
elm["atteints"] = {"n_distincts": int(len(XV)), "essais": list(BANC["gains_elm"]) + list(BANC["essais_hors_J"]["gains_elm"]),
                   "x_min": XV.min(axis=0).tolist(), "x_max": XV.max(axis=0).tolist(),
                   "g_max": float(gv.max()), "pire_MP": float(MV["MP"].min()), "pire_GM": float(MV["GM"].min()),
                   "pire_GB": float(MV["GB"].min()),
                   "pire_module": float(modv.max()), "fc_max": float(MV["FC"].max())}
ORIGINE = {}
for code, liste in list(BANC["gains_elm"].items()) + list(BANC["essais_hors_J"]["gains_elm"].items()):
    for k in liste:
        ORIGINE.setdefault(tuple(np.round(np.array(k) / K_ZN, 12)), set()).add(code)
# charge legere (complementaire) : gains atteints et ensemble (grille 2^(1/16))
leg = []
for c in COINS_LEGERS:
    mA = modules_max(A16, *c)
    mV = modules_max(XV, *c)
    leg.append({"R": c[0], "Vin": c[1], "part_ensemble_stable": float(np.mean(mA < 1.0)),
                "atteints_instables": int(np.sum(mV >= 1.0)),
                "atteints_instables_x": XV[mV >= 1.0].tolist(),
                "atteints_instables_essais": sorted(set().union(*[ORIGINE[tuple(x)] for x in XV[mV >= 1.0]])) if np.any(mV >= 1.0) else [],
                "ZN_module": [z for z in zn if (z["R"], z["Vin"]) == c][0]["module_max"]})
elm["charge_legere"] = leg

# (vi) retard d'un pas Tc sur l'ensemble (grille 2^(1/16))
modr, _ = pire_modules(A16, retard=True)
MR = pire_marges(A16, TERMES_R, TERMES_NY_R)
elm["retard"] = {"pire_module": float(modr.max()), "n_instables": int(np.sum(modr >= 1.0)),
                 "x_pire_module": A16[int(np.argmax(modr))].tolist(),
                 "pire_MP": float(MR["MP"].min()), "x_pire_MP": A16[int(np.argmin(MR["MP"]))].tolist(),
                 "MP_au_pire_MP_ref": float(MR["MP"][i_mp]), "pire_GM": float(MR["GM"].min()),
                 "perte_MP_max": float(np.max(M16["MP"] - MR["MP"])), "perte_MP_min": float(np.min(M16["MP"] - MR["MP"]))}
RES["elm_pid"] = elm
RES["M0_json"] = ENS_JSON["marge_depart_deg"]


# %% ETAPE 5 : figures (sans titre, axes en francais, PDF vectoriel)

# Texte et formules composes par LaTeX avec les polices du memoire (Times : newtxtext et newtxmath) ;
# la virgule decimale des formules vient d'icomma. Prerequis : latex et les paquets newtx et icomma.
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 9, "pdf.fonttype": 42,
                     "text.usetex": True,
                     "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{newtxtext,newtxmath}"
                                            r"\usepackage{icomma}\DeclareUnicodeCharacter{03A9}{\ensuremath{\Omega}}",
                     "font.family": "serif", "axes.grid": True, "grid.alpha": 0.3})


def vf(txt):
    """Graduation au format francais : virgule decimale, vrai signe moins."""
    return txt.replace(".", ",").replace("-", "\u2212")


FR = FuncFormatter(lambda v, _: vf(f"{v:g}"))

# Figure 1 : Bode de la boucle ouverte ZN au nominal
Wb = np.logspace(np.log10(2 * np.pi * 10.0), np.log10(0.98 * np.pi / TC), 6000)
Lb = np.ones(3) @ termes(*NOMINAL, z=np.exp(1j * Wb * TC))
fb = Wb / (2 * np.pi)
mag = 20 * np.log10(np.abs(Lb))
phd = np.degrees(np.unwrap(np.angle(Lb)))
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 5.0), sharex=True)
a1.semilogx(fb, mag, color="#1f4e79", lw=1.4)
a1.axhline(0, color="k", lw=0.7)
a1.axvline(znn["coupure_Hz"], color="#a0a0a0", ls="--", lw=0.8)
a1.axvline(znn["f_180_Hz"], color="#a0a0a0", ls=":", lw=0.8)
mag180 = -znn["marge_gain_dB"]
a1.annotate("", xy=(znn["f_180_Hz"], mag180), xytext=(znn["f_180_Hz"], 0),
            arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.1))
a1.text(znn["f_180_Hz"] * 0.8, mag180 * 0.75, f"MG = {virg(znn['marge_gain_dB'], 1)} dB", color="#c0392b", ha="right", va="center")
for a in (a1, a2):
    a.axvspan(FSW / 2, fb[-1], color="#bbbbbb", alpha=0.25, lw=0)
a1.text(FSW / 2 * 1.08, 40, "f > fs/2", fontsize=8, color="#555555")
a1.set_ylabel("Gain (dB)")
a2.semilogx(fb, phd, color="#1f4e79", lw=1.4)
a2.axhline(-180, color="k", lw=0.7)
a2.axvline(znn["coupure_Hz"], color="#a0a0a0", ls="--", lw=0.8)
a2.axvline(znn["f_180_Hz"], color="#a0a0a0", ls=":", lw=0.8)
a2.annotate("", xy=(znn["coupure_Hz"], -180), xytext=(znn["coupure_Hz"], -180 + znn["marge_phase_deg"]),
            arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.1))
a2.text(znn["coupure_Hz"] * 1.12, -180 + znn["marge_phase_deg"] / 2, f"MP = {virg(znn['marge_phase_deg'], 1)}°",
        color="#c0392b", va="center")
a2.set_ylabel("Phase (°)")
a2.set_xlabel("Fréquence (Hz)")
a2.set_xlim(fb[0], fb[-1])
for a in (a1, a2):
    a.yaxis.set_major_formatter(FR)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "bode_boucle_ouverte_zn_nominal.pdf"))
plt.close(fig)

# Figure 2 : carte des poles ZN (9 coins), vue d'ensemble et agrandissement
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 3.3))
th = np.linspace(0, 2 * np.pi, 600)
cols = {4.0: "#1f77b4", 5.0: "#2ca02c", 6.0: "#d62728"}
marq = {160.0: "o", 200.0: "s", 240.0: "^"}
for a in (a1, a2):
    a.plot(np.cos(th), np.sin(th), color="k", lw=0.8)
    for z in zn[:9]:
        pz = np.array([complex(*p) for p in z["poles"]])
        a.scatter(pz.real, pz.imag, s=22, facecolors="none", edgecolors=cols[z["R"]], marker=marq[z["Vin"]], lw=0.9)
    a.set_xlabel("Partie réelle")
    a.xaxis.set_major_formatter(FR)
    a.yaxis.set_major_formatter(FR)
a1.set_aspect("equal")
a1.set_ylabel("Partie imaginaire")
a2.set_xlim(0.975, 1.003)
a2.set_ylim(-0.03, 0.03)
a2.xaxis.set_major_formatter(FuncFormatter(lambda v_, _: vf(f"{v_:.3f}")))
a2.yaxis.set_major_formatter(FuncFormatter(lambda v_, _: vf(f"{v_:.2f}")))
a1.add_patch(plt.Rectangle((0.975, -0.03), 0.028, 0.06, fill=False, ec="#888888", lw=0.6))
for R, cc in cols.items():
    a1.scatter([], [], facecolors="none", edgecolors=cc, marker="o", label=f"R = {virg(R, 0)} Ω")
for V, mk in marq.items():
    a1.scatter([], [], facecolors="none", edgecolors="k", marker=mk, label=f"Vin = {virg(V, 0)} V")
a1.legend(loc="lower left", fontsize=7, framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "poles_boucle_fermee_zn.pdf"))
plt.close(fig)

# Figure 3 : carte (xP, xI) de la pire marge de phase sur les 9 coins, coupes xD = 1/2, 1, 2
e16 = np.arange(-32, 33) / 16.0
MPg = np.full(len(X16), np.nan)
MPg[np.where(g16 <= 0)[0]] = M16["MP"]
MPg = MPg.reshape(65, 65, 65)
fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.7), sharey=True)
niveaux = np.linspace(M0, np.nanmax(MPg), 12)
for a, kD in zip(axs, (-1.0, 0.0, 1.0)):
    k = int(np.argmin(np.abs(e16 - kD)))
    cs = a.contourf(2.0 ** e16, 2.0 ** e16, MPg[:, :, k].T, levels=niveaux, cmap="viridis")
    a.contour(2.0 ** e16, 2.0 ** e16, g16.reshape(65, 65, 65)[:, :, k].T, levels=[0.0], colors="k", linewidths=0.8)
    if kD == 0.0:
        a.plot(1, 1, "+", color="#c0392b", ms=9, mew=1.6)
    a.set_xscale("log", base=2)
    a.set_yscale("log", base=2)
    a.set_xlabel(r"$\theta_P$")   # multiplicateurs notes theta au chapitre 4
    a.text(0.05, 0.95, rf"$\theta_D$ = {virg(2.0 ** kD, 1)}", transform=a.transAxes, va="top", fontsize=8, color="k")
    a.set_xticks([0.25, 1, 4])
    a.set_yticks([0.25, 0.5, 1, 2, 4])
    a.xaxis.set_major_formatter(FR)
    a.yaxis.set_major_formatter(FR)
    a.xaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    a.yaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    a.grid(False)
axs[0].set_ylabel(r"$\theta_I$")
cb = fig.colorbar(cs, ax=axs, shrink=0.9, pad=0.02)
cb.set_label("Pire marge de phase (°)")
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
fig.savefig(os.path.join(FIG, "elm_marge_phase_ensemble.pdf"), bbox_inches="tight")
plt.close(fig)

# Figure 4 : distributions sur l'ensemble (grille 2^(1/16)) : pire module, pire marge de gain
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 2.8))
a1.hist(mod16, bins=60, color="#1f4e79")
a1.axvline(max(z["module_max"] for z in zn[:9]), color="#c0392b", ls="--", lw=1.0, label="ZN")
a1.set_xlabel("Pire module des pôles (9 coins)")
a1.set_ylabel("Nombre de réglages")
a1.legend()
a2.hist(M16["GM"], bins=60, color="#1f4e79")
a2.axvline(min(z["marge_gain_dB"] for z in zn[:9]), color="#c0392b", ls="--", lw=1.0, label="ZN")
a2.set_xlabel("Pire marge de gain (dB, 9 coins)")
a2.legend()
for a in (a1, a2):
    a.xaxis.set_major_formatter(FR)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "elm_distributions_ensemble.pdf"))
plt.close(fig)


# %% ETAPE 6 : ecriture
with open(os.path.join(ICI, "resultats_stabilite.json"), "w", encoding="utf-8") as f:
    json.dump(RES, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
print(f"Fini en {time.time() - T0:.0f} s.")


# %% ETAPE 7 : resultats_stabilite.md (tableaux)

def v(x, nd=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "sans objet"
    return virg(x, nd)


ZN9 = zn[:9]
ZNL = zn[9:]
E = elm
fm, fp, tb, at, rt, fr_ = E["fine_marges"], E["fine_poles"], E["table"], E["atteints"], E["retard"], E["frontiere"]
cv = RES["convertisseur"]["coins"]
pn = RES["convertisseur"]["poles_nominal"]
L_ = []
w = L_.append
w("# Stabilité linéaire de la boucle à gains figés")
w("")
w("Produit par `stabilite.py`.")
w("")
w("## 1. Modèle")
w("")
w("- Convertisseur : modèle moyen sans pertes, conduction continue, état [iL ; vC], entrée d, sortie vC.")
w("- Discrétisation : bloqueur d'ordre zéro à Tc = 1/220 000 s. Justification (banc `banc_commun.py`, "
  "`simuler`) : d(k) est calculé à l'instant k·Tc à partir de v(k·Tc) et tenu pendant les 120 pas qui suivent.")
w("- Retard de calcul : aucun dans le modèle de référence (le banc applique d(k) à l'instant de la mesure). "
  "Variante de sensibilité : un pas Tc (z⁻¹).")
w("- MLI : gain unité (d comparé à une dent de scie normalisée). Bloc PID commun (éq. 3.4), N = 64 122,9 rad/s, "
  "saturation et anti-emballement ignorés.")
w("- Points : 9 coins R = 4, 5, 6 Ω × Vin = 160, 200, 240 V ; complément hors critère : 25 et 98 Ω.")
w("")
w("## 2. Convertisseur seul")
w("")
w("A = [[0, −1/L], [1/C, −1/(RC)]], B = [Vin/L ; 0], C = [0, 1], D = 0. "
  f"À 5 Ω, 200 V : A = [[0 ; {v(-1 / L_BOB, 0)}], [{v(1 / C_CONV, 1)} ; {v(-1 / (5 * C_CONV), 1)}]], "
  f"B = [{v(200 / L_BOB, 0)} ; 0].")
w("")
w(f"Pôles à 5 Ω : **{v(pn[1])}** et **{v(pn[0])} s⁻¹**. "
  f"ωn = {v(cv[4]['wn'], 1)} rad/s, ζ = {v(cv[4]['zeta'], 3)}.")
w("")
w("| R (Ω) | ζ | Pôles (s⁻¹) | Pôles discrets (module) | écart exp(pTc) |")
w("|---|---|---|---|---|")
for R in (4.0, 5.0, 6.0, 25.0, 98.0):
    c = [x for x in cv if x["R"] == R][0]
    ps = " ; ".join(f"{v(p[0])}" + (f" ± {v(abs(p[1]))} j" if abs(p[1]) > 0 else "") for p in c["poles_s"][:(1 if abs(c['poles_s'][0][1]) > 0 else 2)])
    w(f"| {v(R, 0)} | {v(c['zeta'], 3)} | {ps} | {' ; '.join(v(m, 5) for m in sorted(c['poles_z_module']))} | {c['ecart_zoh']:.1e} |".replace(".", ","))
w("")
w("Les pôles ne dépendent pas de Vin (gain seul).")
w("")
w("## 3. Boucle fermée avec les gains de Ziegler-Nichols")
w("")
w("| R (Ω) | Vin (V) | max |z| | Marge de phase | Coupure (Hz) | Marge de gain (dB) | à (Hz) | MG par les pôles (dB) | min |1+L| | MP avec retard Tc |")
w("|---|---|---|---|---|---|---|---|---|---|")
for z in ZN9:
    w(f"| {v(z['R'], 0)} | {v(z['Vin'], 0)} | {v(z['module_max'], 5)} | {v(z['marge_phase_deg'])}° | {v(z['coupure_Hz'], 0)} | "
      f"{v(z['marge_gain_dB'])} | {v(z['f_180_Hz'], 0)} | {v(z['marge_gain_poles_dB'])} | {v(z['marge_module'], 3)} | "
      f"{v(z['retard_marge_phase_deg'])}° |")
w("")
w("La marge de gain est recoupée de façon indépendante : facteur de gain qui amène un pôle sur le cercle unité "
  "(colonne « par les pôles »). Avec le retard d'un pas, la marge de gain ZN tombe entre "
  f"{v(min(z['retard_marge_gain_dB'] for z in ZN9))} et {v(max(z['retard_marge_gain_dB'] for z in ZN9))} dB et "
  f"le module max reste sous {v(max(z['retard_module_max'] for z in ZN9), 5)}.")
w("")
w("Charge légère (hors critère de l'ensemble, pour mémoire) :")
w("")
w("| R (Ω) | Vin (V) | max |z| | Marge de phase | Coupure (Hz) | Marge de gain à la hausse (dB) | Marge de gain à la baisse (dB) |")
w("|---|---|---|---|---|---|---|")
for z in ZNL:
    inst = z["module_max"] >= 1.0
    w(f"| {v(z['R'], 0)} | {v(z['Vin'], 0)} | {v(z['module_max'], 5)} | {v(z['marge_phase_deg'])}° | {v(z['coupure_Hz'], 0)} | "
      f"{'instable' if inst else v(z['marge_gain_dB'])} | {'instable' if inst else v(z['marge_gain_baisse_dB'])} |")
w("")
w("À charge légère, la résonance du filtre (ζ = 0,29 et 0,074, f0 = 232 Hz) fait passer la phase de L par −180° "
  "à une fréquence où |L| > 1 : la boucle est conditionnellement stable, et une baisse du gain de quelques dB la "
  "déstabilise (marge à la baisse vérifiée par les pôles).")
w("")
w("## 4. ELM-PID à gains figés sur l'ensemble admissible")
w("")
w("Ensemble tel que la loi le voit : g interpolé (trilinéaire en log2) dans `ensemble_gains_elm.mat`, g ≤ 0, "
  "x ∈ [1/4 ; 4]³. Pire cas sur les 9 coins pour chaque x.")
w("")
w("| Jeu de x | Nombre | max |z| | Pire marge de phase | Plus haute coupure (Hz) | Pire marge de gain (dB) | min |1+L| |")
w("|---|---|---|---|---|---|---|")
w(f"| Points admissibles de la table (pas 2^(1/8)) | {tb['n']} | {v(tb['pire_module'], 5)} | {v(tb['pire_MP'], 3)}° | "
  f"{v(tb['fc_max'], 1)} | {v(tb['pire_GM'])} | {v(tb['pire_MM'], 3)} |")
w(f"| Grille 2^(1/16), admis | {fm['n_admis']} | {v(float(mod16.max()), 5)} | {v(fm['pire_MP'], 3)}° | "
  f"{v(fm['fc_max'], 1)} | {v(fm['pire_GM'])} | {v(fm['pire_MM'], 3)} |")
w(f"| Grille 2^(1/32), admis (pôles seuls) | {fp['n_admis']} | {v(fp['pire_module'], 5)} | - | - | - | - |")
w(f"| Gains atteints par l'ELM-PID (S1 à S9, S10) | {at['n_distincts']} | {v(at['pire_module'], 5)} | "
  f"{v(at['pire_MP'], 3)}° | {v(at['fc_max'], 1)} | {v(at['pire_GM'])} | - |")
w("")
w(f"- Aucun x admis n'est instable aux 9 coins : {fp['n_instables']} sur {fp['n_admis']} (grille 2^(1/32)). "
  f"Pire module {v(fp['pire_module'], 5)} en x = {fmt_x(fp['x_pire_module'])}, coin {v(fp['coin_pire_module'][0], 0)} Ω / "
  f"{v(fp['coin_pire_module'][1], 0)} V : un pôle lent (petit xP), pas une oscillation.")
w(f"- Pire marge de phase sur la table : {v(tb['pire_MP'], 3)}° au départ (1 ; 1 ; 1), coin "
  f"{v(tb['coin_pire_MP'][0], 0)} Ω / {v(tb['coin_pire_MP'][1], 0)} V, égale à M0 = {v(RES['M0_json'], 3)}° du JSON.")
w(f"- Entre les points de la table, {fm['n_sous_M0']} points de la grille 2^(1/16) admis par l'interpolation "
  f"ont une marge sous M0 ; la pire vaut {v(fm['pire_MP'], 3)}° (écart {v(fm['pire_MP'] - M0, 3)}°), en "
  f"x = {fmt_x(fm['x_pire_MP'])}. Tous sont à la frontière (g interpolé entre "
  f"{min(fm['g_sous_M0']):.1e} et {max(fm['g_sous_M0']):.1e}). Pire marge à l'intérieur strict : "
  f"{v(fr_['pire_MP_interieur'], 2)}°.")
w(f"- Coupure : au plus {v(fm['fc_max'], 2)} Hz (borne 2 200 Hz ; dépassement de "
  f"{v(fm['fc_max'] - FC_MAXI, 2)} Hz, effet de l'interpolation).")
w(f"- Marge de gain à la hausse : au moins {v(fm['pire_GM'])} dB (x = {fmt_x(fm['x_pire_GM'])}), aucun point "
  f"sous 6 dB ; {v(100 * fm['part_GM_sous_ZN'], 1)} % de l'ensemble est sous la plus petite marge de ZN "
  f"({v(min(z['marge_gain_dB'] for z in ZN9))} dB). {fm['n_conditionnels']} réglages (sur {fm['n_admis']}) sont "
  f"conditionnellement stables à au moins un coin (marge à la baisse au plus {v(fm['pire_GB'])} dB, "
  f"x = {fmt_x(fm['x_pire_GB'])}).")
w(f"- Gains atteints : x entre {fmt_x(at['x_min'])} et {fmt_x(at['x_max'])}, tous admissibles "
  f"(g max {at['g_max']:.1e}), marge de phase au moins {v(at['pire_MP'], 2)}°, marge de gain au moins "
  f"{v(at['pire_GM'])} dB.")
w("")
w("Variante avec un retard d'un pas Tc (sensibilité) : "
  f"aucun x admis instable ({rt['n_instables']}), pire module {v(rt['pire_module'], 5)} ; pire marge de phase "
  f"{v(rt['pire_MP'])}° (x = {fmt_x(rt['x_pire_MP'])}), perte de {v(rt['perte_MP_min'])} à {v(rt['perte_MP_max'])}° "
  f"selon x ; pire marge de gain {v(rt['pire_GM'])} dB.")
w("")
w("Charge légère (hors critère ; l'ensemble ne promet rien ici), grille 2^(1/16) :")
w("")
w("| R (Ω) | Vin (V) | Part de l'ensemble stable | Gains atteints instables | max |z| ZN |")
w("|---|---|---|---|---|")
for c in E["charge_legere"]:
    xs = ", ".join(fmt_x(x) for x in c["atteints_instables_x"]) if c["atteints_instables"] else "aucun"
    w(f"| {v(c['R'], 0)} | {v(c['Vin'], 0)} | {v(100 * c['part_ensemble_stable'], 1)} % | {xs} | {v(c['ZN_module'], 5)} |")
w("")
w("À 98 Ω sous 160 V, ZN lui-même est instable ; les deux autres réglages instables, voisins de ZN, ont été "
  "atteints sur S8b (98 Ω). Le calcul à gains figés ne dit pas si l'essai les a tenus pendant une phase à 160 V.")
w("")
w("## 5. Limites : ce qui est prouvé et ce qui ne l'est pas")
w("")
w("Prouvé, sous les hypothèses du modèle (moyen, linéaire, conduction continue, sans retard de calcul) :")
w("")
w("- la boucle ZN est stable aux 9 coins nominaux, avec 24,9 à 39,5° de marge de phase ;")
w("- tout réglage figé de l'ensemble admissible (table interpolée) est stable aux 9 coins, avec au moins "
  f"{v(fm['pire_MP'], 2)}° de marge de phase et {v(fm['pire_GM'])} dB de marge de gain ; il le reste avec un retard "
  "de calcul d'un pas.")
w("")
w("Non prouvé :")
w("")
w("- la stabilité du système à gains variables : les gains changent toutes les 0,5 ms ; un système linéaire dont "
  "chaque gel est stable peut être instable si ses paramètres varient assez vite. Le calcul ne borne pas la vitesse "
  "d'adaptation (voir le chapitre 4 du mémoire) ;")
w("- la saturation du rapport cyclique, l'anti-emballement, la conduction discontinue, l'échantillonnage de "
  "l'ondulation de découpage (10 échantillons par période) et la modulation multiple dans une période ;")
w(f"- la stabilité à charge légère : à 98 Ω sous 160 V, ZN et {v(100 * (1 - E['charge_legere'][3]['part_ensemble_stable']), 1)} % de l'ensemble sont instables ; l'ensemble "
  "admissible ne contient aucune garantie hors de 4 à 6 Ω ;")
w("- les points entre les coins (R et Vin intermédiaires) et les tolérances de L et C ;")
w("- toute valeur au-delà d'environ fs/2 (11 kHz) : la marge de gain, mesurée vers 26 kHz, en fait partie.")
w("")
w("## 6. Figures (PDF vectoriels, `figures/`)")
w("")
w("- `bode_boucle_ouverte_zn_nominal.pdf` : Bode de L(z) ZN à 5 Ω, 200 V, marges indiquées ; zone grisée au-delà "
  "de fs/2, hors du domaine de validité du modèle moyen.")
w("- `poles_boucle_fermee_zn.pdf` : pôles de la boucle fermée ZN aux 9 coins (disque unité, et agrandissement "
  "autour de z = 1 ; le pôle réel vers 0,72 est proche du pôle du filtre de dérivée, z = 1 − N·Tc = 0,709).")
w("- `elm_marge_phase_ensemble.pdf` : pire marge de phase sur les 9 coins dans le plan (xP, xI), coupes xD = 0,5 ; 1 ; 2, "
  "frontière g = 0 en noir, départ (1 ; 1 ; 1) en croix sur la coupe xD = 1 ; en blanc, réglages non admissibles.")
w("- `elm_distributions_ensemble.pdf` : distributions du pire module des pôles et de la pire marge de gain sur "
  "l'ensemble (grille 2^(1/16)), repère ZN.")
w("")
with open(os.path.join(ICI, "resultats_stabilite.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L_) + "\n")
print("resultats_stabilite.md ecrit.")
