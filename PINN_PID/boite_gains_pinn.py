# =============================================================================
# boite_gains_pinn.py
#
# OBJECTIF
#   Construire la boite des gains admissibles du PINN-PID : le pave dans
#   lequel Adam a le droit de deplacer les gains P, I, D du bloc PID.
#   Gains figes en un point de la boite, la boucle garde, sur la plage
#   nominale, une marge de phase au moins egale a celle du regulateur de
#   depart et une coupure sous fs/10. C'est verifie sur deux grilles (pas
#   2^(1/4) puis 2^(1/8)), avec le modele moyen sans pertes.
#
# LA LOI DONT ON BORNE LES GAINS
#   Bloc "PID Controller" de la base commune : forme Parallel, Forward Euler
#   a Tc pour l'integrateur et le filtre de derivee (pole N = 64 122.9 rad/s) :
#     C(z) = P + I Tc / (z - 1) + D N (z - 1) / (z - 1 + Tc N).
#   Gains de depart : ceux de Ziegler-Nichols du bloc, [P, I, D] =
#   [0.093910, 301.089, 7.3227e-6]. Un point de la boite s'ecrit par ses
#   multiplicateurs (a, b, c) de ces trois gains.
#
# LA REGLE
#   Coins nominaux : R = 4, 5 et 6 ohms ; Vin = 160, 200 et 240 V. Modele
#   moyen sans pertes, Vout/d = Vin / (L C s^2 + (L/R) s + 1), discretise
#   avec un bloqueur d'ordre zero a Tc. Un reglage (a, b, c) est admissible
#   si, aux 9 coins nominaux :
#     1. sa plus petite marge de phase est au moins celle du point de
#        depart (Ziegler-Nichols dans ce bloc : la valeur est imprimee) ;
#     2. sa plus haute frequence de coupure est au plus fs/10 = 2.2 kHz.
#   La boite est le pave qui contient le depart (1, 1, 1), dont tous les
#   points de la grille 2^(k/4) (deux octaves de part et d'autre) sont
#   admissibles, et dont le volume en echelle logarithmique est le plus
#   grand.
#
# CONTROLES
#   - la boite est reverifiee sur une grille deux fois plus fine (2^(k/8)) ;
#   - marges a charge legere (25 et 98 ohms) aux coins de la boite, pour
#     information.
#
# FICHIERS
#   Lus : banc_commun.py, scenarios_communs.json et scenario_S*.mat (lus par
#     le banc commun ; aucun essai n'est simule ici).
#   Ecrit : boite_gains_pinn.json.
#
# COMMENT LANCER CE SCRIPT
#   python boite_gains_pinn.py      (moins d'une minute)
#
# ORDRE D'EXECUTION
#   1. ce script ; 2. choix_correction_pinn.py ; 3. banc_pinn_pid.py.
# =============================================================================


# %% ETAPE 0 : bibliotheques et definitions du banc commun

import os                                     # chemins de fichiers
import json                                   # fichier de resultats
import time                                   # duree d'execution
import numpy as np                            # calcul numerique
from scipy.signal import cont2discrete        # discretisation du modele moyen

try:                                          # dossier du script
    DOSSIER_PINN = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PINN = os.getcwd()

# Le banc commun est execute jusqu'a sa marque de fin des definitions :
# constantes du circuit (L_BOB, C_CONV, FSW), Tc et PID de reference.
with open(os.path.join(DOSSIER_PINN, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py n'a pas la marque de fin des definitions.")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])

T_DEBUT = time.time()                         # pour la duree totale
K_ZN = np.array([PID_P, PID_I, PID_D])       # gains de depart : Ziegler-Nichols du bloc PID [P, I, D]
FC_MAXI = FSW / 10.0                          # critere 2 : frequence de coupure maximale (Hz)
PAS = 2.0 ** 0.25                             # pas de la grille (multiplicatif)
OCTAVES = 2                                   # etendue de la grille de part et d'autre du depart
COINS_NOMINAUX = [(R, V) for R in (4.0, 5.0, 6.0) for V in (160.0, 200.0, 240.0)]
COINS_LEGERS = [(R, V) for R in (25.0, 98.0) for V in (160.0, 200.0, 240.0)]
print(f"Gains de depart (Ziegler-Nichols, bloc PID parallele) : P = {K_ZN[0]:.6f}, I = {K_ZN[1]:.6e}, "
      f"D = {K_ZN[2]:.6e}")


# %% ETAPE 1 : marge de phase et frequence de coupure d'un reglage

W = np.logspace(np.log10(2 * np.pi * 20.0), np.log10(np.pi / TC * 0.999), 8000)   # pulsations (rad/s)
Z = np.exp(1j * W * TC)                       # points du cercle unite


def reponse_buck(R, V):
    """Reponse frequentielle du modele moyen du Buck discretise (bloqueur
    d'ordre zero a Tc) au coin (R, V)."""
    num_d, den_d, _ = cont2discrete(([V], [L_BOB * C_CONV, L_BOB / R, 1.0]), TC, method="zoh")
    return np.polyval(np.squeeze(num_d), Z) / np.polyval(den_d, Z)


def preparer(coins):
    """Les trois termes du gain de boucle pour des gains unitaires (multiplicateur 1),
    a chaque coin : terme proportionnel, integral et derive."""
    G = np.array([reponse_buck(R, V) for R, V in coins])          # une ligne par coin
    terme_p = K_ZN[0] * G
    terme_i = K_ZN[1] * TC / (Z - 1.0) * G                         # bloc parallele : I Tc / (z - 1)
    terme_d = K_ZN[2] * PID_N * (Z - 1.0) / (Z - 1.0 + TC * PID_N) * G
    return terme_p, terme_i, terme_d


TERMES_NOM = preparer(COINS_NOMINAUX)
TERMES_LEG = preparer(COINS_LEGERS)


def frequentiel(a, b, c, termes):
    """Plus petite marge de phase (degres) et plus haute frequence de
    coupure (Hz) sur les coins de 'termes', pour le reglage (a, b, c). Le
    passage du gain par 1 est interpole lineairement en echelle log."""
    tp, ti, td = termes
    boucle = a * tp + b * ti + c * td                             # gain de boucle a chaque coin
    lm = np.log(np.abs(boucle))                                   # log du module
    ph = np.unwrap(np.angle(boucle), axis=1)                      # phase deroulee
    lignes, j = np.where(np.diff(np.sign(lm), axis=1) != 0)       # le module traverse 1
    if len(j) == 0:
        return 180.0, 0.0
    f = lm[lignes, j] / (lm[lignes, j] - lm[lignes, j + 1])       # position du passage entre j et j+1
    phase = ph[lignes, j] + f * (ph[lignes, j + 1] - ph[lignes, j])
    wc = W[j] * (W[j + 1] / W[j]) ** f
    marges = (np.degrees(phase) + 360.0) % 360.0 - 180.0
    return float(marges.min()), float(wc.max() / (2.0 * np.pi))


def marge_de_gain(a, b, c, termes):
    """Plus petite marge de gain (dB) sur les coins de 'termes' : distance
    du module a 1 la ou la phase traverse -180 degres (a 360 degres pres).
    Pour information seulement : elle n'entre pas dans la regle."""
    tp, ti, td = termes
    boucle = a * tp + b * ti + c * td
    phase = np.degrees(np.unwrap(np.angle(boucle), axis=1))
    tours = np.floor((phase + 180.0) / 360.0)                     # change quand la phase traverse -180
    lignes, j = np.where(np.diff(tours, axis=1) != 0)
    if len(j) == 0:
        return float("inf")
    return float(np.min(-20.0 * np.log10(np.abs(boucle[lignes, j]))))


MARGE_DEPART, FC_DEPART = frequentiel(1.0, 1.0, 1.0, TERMES_NOM)
MARGE_DEPART_LEG, _ = frequentiel(1.0, 1.0, 1.0, TERMES_LEG)
print(f"Point de depart : marge minimale {MARGE_DEPART:.2f} deg et coupure maximale {FC_DEPART:.0f} Hz aux coins "
      f"nominaux ; marge minimale {MARGE_DEPART_LEG:.1f} deg a charge legere.")


def admissible(a, b, c):
    """Criteres 1 et 2 aux 9 coins nominaux."""
    marge, fc = frequentiel(a, b, c, TERMES_NOM)
    return marge >= MARGE_DEPART - 1e-6 and fc <= FC_MAXI


# %% ETAPE 2 : grille et recherche de la plus grande boite

titre = lambda txt: print("\n" + "=" * 78 + "\n " + txt + "\n" + "=" * 78)
titre("ETAPE 2 : grille des reglages admissibles et plus grande boite")
EXPOSANTS = np.arange(-4 * OCTAVES, 4 * OCTAVES + 1)              # multiplicateurs PAS^k
MULT = PAS ** EXPOSANTS
n = len(MULT)
i0 = int(np.where(EXPOSANTS == 0)[0][0])      # indice du multiplicateur 1
OK = np.zeros((n, n, n), bool)                # OK[i, j, k] : reglage (MULT[i], MULT[j], MULT[k]) admissible
for i, a in enumerate(MULT):
    for j, b in enumerate(MULT):
        for k, c in enumerate(MULT):
            OK[i, j, k] = admissible(a, b, c)
print(f"  {n ** 3} reglages, {int(OK.sum())} admissibles ({time.time() - T_DEBUT:.0f} s).")

# Sommes cumulees des reglages NON admissibles : le nombre de points
# non admissibles d'un pave s'obtient en 8 lectures.
S = np.zeros((n + 1, n + 1, n + 1), int)
S[1:, 1:, 1:] = (~OK).astype(int).cumsum(0).cumsum(1).cumsum(2)


def non_admissibles(i1, i2, j1, j2, k1, k2):
    """Nombre de points non admissibles dans le pave d'indices [i1, i2] x [j1, j2] x [k1, k2] (bornes comprises)."""
    return (S[i2 + 1, j2 + 1, k2 + 1] - S[i1, j2 + 1, k2 + 1] - S[i2 + 1, j1, k2 + 1] - S[i2 + 1, j2 + 1, k1]
            + S[i1, j1, k2 + 1] + S[i1, j2 + 1, k1] + S[i2 + 1, j1, k1] - S[i1, j1, k1])


paves = []                                    # (volume en pas de grille, indices) des paves admissibles
for i1 in range(i0 + 1):
    for i2 in range(i0, n):
        for j1 in range(i0 + 1):
            for j2 in range(i0, n):
                for k1 in range(i0 + 1):
                    for k2 in range(i0, n):
                        if non_admissibles(i1, i2, j1, j2, k1, k2) == 0:
                            paves.append(((i2 - i1) * (j2 - j1) * (k2 - k1), (i1, i2, j1, j2, k1, k2)))
paves.sort(key=lambda p: -p[0])
if not paves or paves[0][0] == 0:
    raise RuntimeError("Aucune boite de volume non nul ne contient le point de depart.")
VOLUME_MAX = paves[0][0]
meilleurs = [p for p in paves if p[0] == VOLUME_MAX]
print(f"  {len(paves)} paves admissibles contiennent le point de depart. Volume maximal : {VOLUME_MAX} "
      f"(en pas de grille au cube), atteint par {len(meilleurs)} pave(s).")
for v, (i1, i2, j1, j2, k1, k2) in paves[:5]:
    print(f"    volume {v:4d} : a de {MULT[i1]:.3f} a {MULT[i2]:.3f}, b de {MULT[j1]:.3f} a {MULT[j2]:.3f}, "
          f"c de {MULT[k1]:.3f} a {MULT[k2]:.3f}")
if len(meilleurs) > 1:
    raise RuntimeError("Plusieurs boites de volume maximal : la regle ne suffit pas a choisir.")
i1, i2, j1, j2, k1, k2 = meilleurs[0][1]
BORNES_MULT = np.array([[MULT[i1], MULT[j1], MULT[k1]], [MULT[i2], MULT[j2], MULT[k2]]])   # [min ; max]
K_MIN, K_MAX = K_ZN * BORNES_MULT[0], K_ZN * BORNES_MULT[1]
print(f"\n  Boite retenue (multiplicateurs) : a de {BORNES_MULT[0, 0]:.4f} a {BORNES_MULT[1, 0]:.4f}, "
      f"b de {BORNES_MULT[0, 1]:.4f} a {BORNES_MULT[1, 1]:.4f}, c de {BORNES_MULT[0, 2]:.4f} a {BORNES_MULT[1, 2]:.4f}")
print(f"  En gains : P de {K_MIN[0]:.6f} a {K_MAX[0]:.6f} ; I de {K_MIN[1]:.4f} a {K_MAX[1]:.4f} ; "
      f"D de {K_MIN[2]:.6e} a {K_MAX[2]:.6e}")
borne_grille = [bool(x) for x in (i1 == 0, i2 == n - 1, j1 == 0, j2 == n - 1, k1 == 0, k2 == n - 1)]
print("  Bornes fixees par la limite de deux octaves : "
      + (", ".join(nom for nom, x in zip(["a min", "a max", "b min", "b max", "c min", "c max"], borne_grille) if x)
         or "aucune"))


# %% ETAPE 3 : controles de la boite

titre("ETAPE 3 : controles de la boite")
# 1. Grille deux fois plus fine a l'interieur de la boite
fin = [PAS ** (np.arange(EXPOSANTS[lo] * 2, EXPOSANTS[hi] * 2 + 1) / 2.0)
       for lo, hi in ((i1, i2), (j1, j2), (k1, k2))]
pire_marge, pire_fc = 180.0, 0.0
for a in fin[0]:
    for b in fin[1]:
        for c in fin[2]:
            m_, f_ = frequentiel(a, b, c, TERMES_NOM)
            pire_marge, pire_fc = min(pire_marge, m_), max(pire_fc, f_)
n_fin = len(fin[0]) * len(fin[1]) * len(fin[2])
print(f"  Grille fine (pas 2^(1/8), {n_fin} reglages) : marge minimale {pire_marge:.2f} deg (seuil "
      f"{MARGE_DEPART:.2f}), coupure maximale {pire_fc:.0f} Hz (seuil {FC_MAXI:.0f}).")
if pire_marge < MARGE_DEPART - 0.5 or pire_fc > FC_MAXI * 1.02:
    raise RuntimeError("La boite viole les criteres entre les points de la grille.")

# 2. Les huit coins de la boite, aux coins nominaux et a charge legere
coins_boite = []
print("  Coins de la boite (multiplicateurs a, b, c) : marge de phase, coupure et marge de gain nominales ; marge de")
print("  phase a charge legere")
for a in BORNES_MULT[:, 0]:
    for b in BORNES_MULT[:, 1]:
        for c in BORNES_MULT[:, 2]:
            m_n, f_n = frequentiel(a, b, c, TERMES_NOM)
            m_l, f_l = frequentiel(a, b, c, TERMES_LEG)
            mg_n = marge_de_gain(a, b, c, TERMES_NOM)
            coins_boite.append({"a": float(a), "b": float(b), "c": float(c), "marge_nominale_deg": m_n,
                                "coupure_nominale_Hz": f_n, "marge_gain_nominale_dB": mg_n,
                                "marge_legere_deg": m_l, "coupure_legere_Hz": f_l})
            print(f"    ({a:5.3f}, {b:5.3f}, {c:5.3f}) : {m_n:5.1f} deg, {f_n:5.0f} Hz, {mg_n:5.1f} dB ; "
                  f"{m_l:6.1f} deg a charge legere")
SENS = []                                     # ce que la boite permet a chaque gain, depuis le depart
for nom, lo, hi in zip(("P", "I", "D"), BORNES_MULT[0], BORNES_MULT[1]):
    SENS.append(f"{nom} " + ("fixe" if lo == hi == 1.0 else "ne peut que monter" if lo == 1.0 else
                             "ne peut que baisser" if hi == 1.0 else "peut monter ou baisser"))
print("  Lecture : depuis le depart (1, 1, 1), " + ", ".join(SENS) + ".")


# %% ETAPE 4 : ecriture du fichier

resultat = {"loi": "bloc PID parallele de la base commune (Forward Euler, derivee filtree N), C(z) = P + I Tc/(z-1) + D N (z-1)/(z-1+Tc N)",
            "K_depart": K_ZN.tolist(), "K_min": K_MIN.tolist(), "K_max": K_MAX.tolist(),
            "multiplicateurs_min": BORNES_MULT[0].tolist(), "multiplicateurs_max": BORNES_MULT[1].tolist(),
            "critere": {"coins_nominaux": COINS_NOMINAUX, "marge_minimale_deg": MARGE_DEPART,
                        "coupure_maximale_Hz": FC_MAXI, "pas_grille": PAS, "octaves": OCTAVES},
            "depart": {"marge_nominale_deg": MARGE_DEPART, "coupure_nominale_Hz": FC_DEPART,
                       "marge_legere_deg": MARGE_DEPART_LEG},
            "grille_fine": {"marge_minimale_deg": pire_marge, "coupure_maximale_Hz": pire_fc, "reglages": n_fin},
            "coins_boite": coins_boite, "bornes_fixees_par_la_grille": borne_grille,
            "volume_max": int(VOLUME_MAX), "paves_admissibles": len(paves), "sens_depuis_depart": SENS}
with open(os.path.join(DOSSIER_PINN, "boite_gains_pinn.json"), "w", encoding="utf-8") as f:
    json.dump(resultat, f, indent=1)
print(f"\n  boite_gains_pinn.json ecrit. Duree totale : {time.time() - T_DEBUT:.0f} s.")
