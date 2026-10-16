# =============================================================================
# banc_fuzzy_pid.py
#
# OBJECTIF
#   Faire tourner le Fuzzy-PID sur le banc commun, sur les onze essais, a cote
#   de Ziegler-Nichols (meme PID, gains fixes) et du meilleur PID fige, et
#   ecrire ce que Simulink doit retrouver :
#     - predictions_banc_fuzzy_pid.json : resultats attendus sur les onze
#       essais, au format des autres fichiers de predictions (lu par
#       Simuler_Fuzzy_PID.m et Construction_Fuzzy_PID.m) ;
#     - fuzzy_pid_reglages.mat : les nombres de l'ordonnanceur (plages,
#       echelles, tables), lus par ordonnanceur_flou_zhao.m, pour que le bloc
#       et ce banc utilisent les memes nombres ;
#     - reference_rejeu_fuzzy_pid.mat : erreur et gains pas par pas sur
#       quatre essais, et le systeme flou seul sur une grille, que
#       Tester_Fuzzy_PID_Rejeu.m rejoue dans le bloc MATLAB.
#
# LA METHODE (ordonnancement flou de Zhao, Tomizuka et Isaka, 1993)
#   A chaque periode Tc = 1/220 000 s :
#   1. e(k) = consigne - mesure, De(k) = e(k) - e(k-1) (0 au premier pas) ;
#      e' = e / E_M, De' = De / DE_M ; sept ensembles triangulaires par
#      entree (NB ... PB), extremites prolongees.
#   2. 49 regles (tableaux I, II, III de l'article) : poids mu_i = produit
#      des appartenances ; K'p = somme mu_i x_C(mu_i), K'd idem, alpha =
#      somme mu_i alpha_i, avec x_Small(mu) = exp(-4 mu), x_Big(mu) =
#      1 - exp(-4 mu).
#   3. Kp = Kp_min + (Kp_max - Kp_min) K'p ; Kd = Kd_min + (Kd_max - Kd_min) K'd ;
#      Ki = Kp^2 / (alpha Kd). Ce sont les champs P, I, D du bloc PID.
#   4. Le PID est le bloc "PID Controller" de Ziegler-Nichols (forme
#      Parallel, Forward Euler, filtre N, saturation [0.01 ; 0.99], clamping)
#      avec P, I, D recus a ce pas : memes operations que PIDClassique du
#      banc commun.
#
# FICHIERS NECESSAIRES (dans le meme dossier que ce script)
#   banc_commun.py, scenarios_communs.json, scenario_S*.mat,
#   predictions_banc_pid_fige.json.
#
# CE QUE PRODUIT CE SCRIPT
#   predictions_banc_fuzzy_pid.json, fuzzy_pid_reglages.mat,
#   reference_rejeu_fuzzy_pid.mat, banc_fuzzy_pid_resultats.json,
#   banc_fuzzy_pid.png (tension, gains et courant sur S1, S8b et S9),
#   fuzzy_pid_surfaces.png (K'p, K'd et alpha en fonction de e' et De').
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT
#   python banc_fuzzy_pid.py      (quelques minutes)
#
# ORDRE D'EXECUTION
#   1. ce script ; 2. Tester_Fuzzy_PID_Rejeu.m ; 3. Construction_Fuzzy_PID.m ;
#   4. Simuler_Fuzzy_PID.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques, reglages et definitions du banc commun

import os                                     # chemins de fichiers
import json                                   # lecture et ecriture au format texte
import math                                   # exponentielle (meme fonction scalaire que MATLAB)
import time                                   # duree d'execution
import numpy as np                            # calcul numerique
from scipy.io import savemat                  # fichiers .mat
import matplotlib.pyplot as plt               # figures

try:                                          # dossier du script
    DOSSIER_FUZZY = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_FUZZY = os.getcwd()

with open(os.path.join(DOSSIER_FUZZY, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC or "iL_max_dem_A" not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py ne contient pas les definitions attendues (marque de fin, observation du courant).")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])  # constantes, essais, Circuit, simuler, grandeurs, PIDClassique

T_DEBUT = time.time()                         # pour la duree totale
ESSAIS = {sc["code"]: sc for sc in SCENARIOS}  # les onze essais, par code
CODES_IAE = ["S1", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]   # IAE de 30 ms a la fin
CODES_DEM = ["S1", "S8a", "S8b"]              # IAE du demarrage (5, 25 et 98 ohms)
K30 = int(round(0.03 / TC))                   # instant 30 ms (numero de pas Tc)
CODES_REJEU = ["S1", "S7a", "S8b", "S9"]      # essais rejoues dans le bloc MATLAB

# Ku et Tu equivalents (regle du point critique appliquee a l'envers aux gains de Ziegler-Nichols)
KU = PID_P / 0.6                              # Kp = 0.6 Ku
TU = 2.0 * PID_P / PID_I                      # Ti = P / I = 0.5 Tu
KP_MIN, KP_MAX = 0.32 * KU, 0.6 * KU          # eq. 12
KD_MIN, KD_MAX = 0.08 * KU * TU, 0.15 * KU * TU
# Echelles de e et De
E_M = 100.0                                   # hauteur de l'echelon de consigne (V)
PENTE_BO = 81937.0                            # pente maximale de la reponse indicielle en boucle ouverte (V/s)
DE_M = PENTE_BO * TC                          # (V par periode Tc)

# Tableaux I, II, III de l'article (lignes e = NB..PB, colonnes De = NB..PB) ; 1 = Big, 0 = Small
B, S = 1, 0
TABLE_KP = np.array([[B, B, B, B, B, B, B],
                     [S, B, B, B, B, B, S],
                     [S, S, B, B, B, S, S],
                     [S, S, S, B, S, S, S],
                     [S, S, B, B, B, S, S],
                     [S, B, B, B, B, B, S],
                     [B, B, B, B, B, B, B]])
TABLE_KD = np.array([[S, S, S, S, S, S, S],
                     [B, B, S, S, S, B, B],
                     [B, B, B, S, B, B, B],
                     [B, B, B, B, B, B, B],
                     [B, B, B, S, B, B, B],
                     [B, B, S, S, S, B, B],
                     [S, S, S, S, S, S, S]])
TABLE_ALPHA = np.array([[2, 2, 2, 2, 2, 2, 2],
                        [3, 3, 2, 2, 2, 3, 3],
                        [4, 3, 3, 2, 3, 3, 4],
                        [5, 4, 3, 3, 3, 4, 5],
                        [4, 3, 3, 2, 3, 3, 4],
                        [3, 3, 2, 2, 2, 3, 3],
                        [2, 2, 2, 2, 2, 2, 2]], dtype=float)
# Controle avec le texte de l'article (points a1 et b1 de la fig. 4)
if not (TABLE_KP[6, 3] == B and TABLE_KD[6, 3] == S and TABLE_ALPHA[6, 3] == 2
        and TABLE_KP[3, 0] == S and TABLE_KD[3, 0] == B and TABLE_ALPHA[3, 0] == 5):
    raise RuntimeError("Les tables ne redonnent pas les deux regles citees dans le texte de l'article.")
for T in (TABLE_KP, TABLE_KD, TABLE_ALPHA):   # tables symetriques par rapport au centre, comme dans l'article
    if not np.array_equal(T, T[::-1, ::-1]):
        raise RuntimeError("Table non symetrique : erreur de recopie.")

K_ZN = np.array([PID_P, PID_I, PID_D])        # Ziegler-Nichols (champs P, I, D du bloc)
print(f"Ku = {KU:.6f}, Tu = {TU * 1e3:.5f} ms ; Kp de {KP_MIN:.6f} a {KP_MAX:.6f} ({KP_MIN / PID_P:.3f} a "
      f"{KP_MAX / PID_P:.3f} fois P) ; Kd de {KD_MIN:.4e} a {KD_MAX:.4e} ({KD_MIN / PID_D:.3f} a {KD_MAX / PID_D:.3f} fois D).")
print(f"Echelles : E_M = {E_M:.1f} V, DE_M = {DE_M:.5f} V par periode Tc.")

# Meilleur PID fige (reference oracle)
with open(os.path.join(DOSSIER_FUZZY, "predictions_banc_pid_fige.json"), encoding="utf-8") as f:
    GAINS_FIGE = json.load(f)["gains"]


def titre(texte):
    """Bandeau de titre dans la console."""
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78)


# %% ETAPE 1 : l'ordonnanceur flou (copie de ordonnanceur_flou_zhao.m) et le PID

def appartenances(x):
    """Degres d'appartenance de x (deja normalise) aux sept ensembles
    NB NM NS ZO PS PM PB : triangles centres en -1, -2/3, ..., 1, de
    demi-largeur 1/3, extremites prolongees (memes formules que le bloc)."""
    xs = min(max(x, -1.0), 1.0)
    return [max(0.0, 1.0 - abs(3.0 * xs - (j - 3))) for j in range(7)]


def zhao(e, de, echelle_e=1.0, echelle_de=1.0):
    """K'p, K'd et alpha de l'article pour (e, De) en volts (eq. 6 a 10).
    Les regles sont parcourues ligne par ligne, celles de poids nul sautees,
    dans le meme ordre que le bloc MATLAB."""
    mE = appartenances(e / (E_M * echelle_e))
    mD = appartenances(de / (DE_M * echelle_de))
    kp = kd = al = 0.0
    for i in range(7):
        if mE[i] > 0.0:
            for j in range(7):
                if mD[j] > 0.0:
                    w = mE[i] * mD[j]                    # eq. 8 (produit)
                    q = math.exp(-4.0 * w)               # x_Small(w) ; x_Big(w) = 1 - q
                    kp += w * (1.0 - q if TABLE_KP[i, j] else q)    # eq. 10a
                    kd += w * (1.0 - q if TABLE_KD[i, j] else q)    # eq. 10b
                    al += w * TABLE_ALPHA[i, j]                     # eq. 10c
    return kp, kd, al


def gains_zhao(e, de, **ech):
    """Champs P, I, D du bloc PID (eq. 11)."""
    kp, kd, al = zhao(e, de, **ech)
    Kp = KP_MIN + (KP_MAX - KP_MIN) * kp
    Kd = KD_MIN + (KD_MAX - KD_MIN) * kd
    return Kp, Kp * Kp / (al * Kd), Kd


class FuzzyPID:
    """Le bloc PID de Ziegler-Nichols (memes operations que PIDClassique du
    banc commun) dont P, I et D sont donnes a chaque pas par l'ordonnanceur
    de Zhao. gains_fixes = True : gains de Ziegler-Nichols (controle)."""

    def __init__(self, gains_fixes=False, **ech):
        self.xI = PID_CI_INTEGRATEUR          # etat de l'integrateur
        self.xF = PID_CI_FILTRE               # etat du filtre de derivee
        self.e1 = None                        # erreur du pas precedent
        self.fixes, self.ech = gains_fixes, ech
        self.K_pas = []                       # gains [P I D] utilises a chaque pas
        self.zh_pas = []                      # [K'p K'd alpha] a chaque pas

    def pas(self, e):
        de = 0.0 if self.e1 is None else e - self.e1     # De(0) = 0
        self.e1 = e
        if self.fixes:
            P, I, D = PID_P, PID_I, PID_D
            self.zh_pas.append((float("nan"),) * 3)
        else:
            kp, kd, al = zhao(e, de, **self.ech)
            P = KP_MIN + (KP_MAX - KP_MIN) * kp
            D = KD_MIN + (KD_MAX - KD_MIN) * kd
            I = P * P / (al * D)
            self.zh_pas.append((kp, kd, al))
        self.K_pas.append((P, I, D))
        derivee = PID_N * (D * e - self.xF)              # terme derive filtre
        u = P * e + self.xI + derivee                    # somme des trois termes
        u_sat = min(max(u, D_MIN), D_MAX)
        entree_I = I * e
        bloque = (u != u_sat) and (np.sign(entree_I) == np.sign(u - u_sat))
        if not bloque:
            self.xI += TC * entree_I
        self.xF += TC * derivee
        return u_sat


# Controle 1 : points du systeme flou calcules a la main
q1 = 1.0 - math.exp(-4.0)                     # x_Big(1)
controles = [((0.0, 0.0), (q1, q1, 3.0)),      # (ZO, ZO) seule : Big, Big, 3
             ((100.0, 0.0), (q1, math.exp(-4.0), 2.0)),         # (PB, ZO) : Big, Small, 2 (point a1)
             ((0.0, -DE_M), (math.exp(-4.0), q1, 5.0)),         # (ZO, NB) : Small, Big, 5 (point b1)
             ((50.0, 0.0), (2 * 0.5 * (1 - math.exp(-2.0)), 2 * 0.5 * math.exp(-2.0), 2.0))]
for (e_, de_), att in controles:          # (50 V = e' 0.5 : moitie PS, moitie PM ; De ZO)
    obt = zhao(e_, de_)
    if max(abs(a - b) for a, b in zip(obt, att)) > 1e-12:
        raise RuntimeError(f"Systeme flou : ({e_}, {de_}) donne {obt}, attendu {att}.")
print("Controle du systeme flou : quatre points calcules a la main retrouves.")
Kp0, Ki0, Kd0 = gains_zhao(0.0, 0.0)
print(f"Gains au repos (ZO, ZO) : Kp {Kp0 / PID_P:.3f} P, Ki {Ki0 / PID_I:.3f} I, Kd {Kd0 / PID_D:.3f} D.")


# %% ETAPE 2 : evaluation d'un regulateur sur les onze essais

def evaluer(fabrique, codes=None, essais=None):
    """Simule un nouveau regulateur (fabrique()) sur chaque essai et renvoie,
    par essai : grandeurs du banc commun, IAE du demarrage, retours, et le
    regulateur (pour ses gains). Essais : les onze (ESSAIS) par defaut."""
    essais = essais or ESSAIS
    res = {}
    for code in (codes or list(essais)):
        reg = fabrique()
        sim = simuler(reg, essais[code])
        o = grandeurs(sim, essais[code])
        e = sim["consigne"] - sim["v"]
        res[code] = {"o": o, "IAE_dem": float(np.sum(np.abs(e[:K30])) * TC),
                     "revenus": bool(all(ev["revenu"] for ev in o["evenements"])), "reg": reg, "sim": sim}
    return res


def termes(res):
    """Les 13 termes du cout du meilleur PID fige."""
    return np.array([res[c]["o"]["IAE"] for c in CODES_IAE] + [res[c]["IAE_dem"] for c in CODES_DEM])


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


titre("ETAPE 2 : les references")
REF = {}
REF["Ziegler-Nichols"] = evaluer(lambda: PIDClassique())
REF["meilleur PID fige"] = evaluer(lambda: PIDParallele(GAINS_FIGE["P"], GAINS_FIGE["I"], GAINS_FIGE["D"]))
T_ZN = termes(REF["Ziegler-Nichols"])
J_FIGE = float(np.mean(termes(REF["meilleur PID fige"]) / T_ZN))
print(f"  Meilleur PID fige : J = {J_FIGE:.3f} (recherche_pid_fige.py : 0.738).")

# Controle 2 : gains forces a Ziegler-Nichols = PID classique du banc a l'identique
fige = evaluer(lambda: FuzzyPID(gains_fixes=True), ["S1", "S8b"])
ecart_zn = max(float(np.max(np.abs(fige[c]["sim"]["d"] - REF["Ziegler-Nichols"][c]["sim"]["d"]))) for c in fige)
print(f"  Controle : Fuzzy-PID a gains forces contre PID classique, ecart maximal sur d : {ecart_zn:.1e}.")
if ecart_zn > 0.0:
    raise RuntimeError("A gains forces, le Fuzzy-PID ne redonne pas le PID classique du banc.")


# %% ETAPE 3 : le Fuzzy-PID sur les onze essais

titre("ETAPE 3 : le Fuzzy-PID (Zhao) sur les onze essais")


def bilan(nom, res):
    """Cout J, retours et une ligne par essai."""
    J = float(np.mean(termes(res) / T_ZN))
    non = [c for c in CODES_IAE if not res[c]["revenus"]]
    print(f"\n  {nom} : J = {J:.3f} ; {'tous les evenements reviennent' if not non else 'pas revenus : ' + ', '.join(non)}")
    for code, r in res.items():
        o = r["o"]
        ligne = (f"    {code:4s} depassement {o['depassement_pct']:5.1f} %, etabli {o['t_etab_ms']:5.2f} ms, "
                 f"IAE {o['IAE'] * 1e3:7.2f} mV.s, efficace {o['e_eff_V'] * 1e3:7.1f} mV, |dd| {o['dd_moyen']:.4f}, "
                 f"iL crete {o['iL_max_dem_A']:5.1f} A")
        if isinstance(r["reg"], FuzzyPID) and not r["reg"].fixes:
            K = np.array(r["reg"].K_pas)[K30:] / K_ZN
            ligne += (f" ; apres 30 ms : Kp {K[:, 0].min():.2f}-{K[:, 0].max():.2f}, Ki {K[:, 1].min():.2f}-"
                      f"{K[:, 1].max():.2f}, Kd {K[:, 2].min():.2f}-{K[:, 2].max():.2f} (x ZN)")
        print(ligne)
    return J


J_ZN = bilan("Ziegler-Nichols (meme PID, gains fixes)", REF["Ziegler-Nichols"])
FZ = evaluer(lambda: FuzzyPID())
J_FZ = bilan("Fuzzy-PID", FZ)

print("\n  Evenements, Fuzzy-PID (Ziegler-Nichols entre parentheses) :")
for code in ("S2", "S3", "S4", "S5", "S6", "S8a", "S8b", "S9"):
    for ev, ev0 in zip(FZ[code]["o"]["evenements"], REF["Ziegler-Nichols"][code]["o"]["evenements"]):
        etat = ("reste dans la bande" if ev["reste_dans_bande"] else
                f"retour en {ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU")
        print(f"    {code:4s} a {ev['t_ms']:6.1f} ms : ecart max {ev['e_max_V']:6.2f} V ({ev0['e_max_V']:6.2f}), "
              f"IAE {ev['IAE'] * 1e3:7.2f} mV.s ({ev0['IAE'] * 1e3:7.2f}), crete a crete en fin "
              f"{ev['e_crete_a_crete_fin_V']:5.2f} V ({ev0['e_crete_a_crete_fin_V']:5.2f}), {etat}")


# %% ETAPE 4 : sensibilite (pas des candidats)

titre("ETAPE 4 : sensibilite (L ou C a +-0.1 % ; echelles E_M et DE_M x0.5 et x2)")
L_NOM, C_NOM = L_BOB, C_CONV                  # valeurs nominales, retablies a la fin
SENSIBILITE = []
print(f"  {'variante':16s} {'J':>6s}  {'non revenus':12s} {'S1 dep. %':>9s} {'S7a eff':>8s} {'S8b IAE':>8s} {'S9 IAE':>8s}")
for nom_p, fl, fc, ech in (("nominal", 1.0, 1.0, {}), ("L - 0.1 %", 0.999, 1.0, {}), ("L + 0.1 %", 1.001, 1.0, {}),
                           ("C - 0.1 %", 1.0, 0.999, {}), ("C + 0.1 %", 1.0, 1.001, {}),
                           ("E_M x0.5", 1.0, 1.0, {"echelle_e": 0.5}), ("E_M x2", 1.0, 1.0, {"echelle_e": 2.0}),
                           ("DE_M x0.5", 1.0, 1.0, {"echelle_de": 0.5}), ("DE_M x2", 1.0, 1.0, {"echelle_de": 2.0})):
    L_BOB, C_CONV = L_NOM * fl, C_NOM * fc   # le banc lit ces constantes a chaque nouveau circuit
    r = FZ if nom_p == "nominal" else evaluer(lambda: FuzzyPID(**ech))
    J_v = float(np.mean(termes(r) / T_ZN))
    non = [c for c in CODES_IAE if not r[c]["revenus"]]
    SENSIBILITE.append({"variante": nom_p, "J": J_v, "non_revenus": non,
                        "S1_depassement_pct": r["S1"]["o"]["depassement_pct"], "S7a_e_eff_V": r["S7a"]["o"]["e_eff_V"],
                        "S8b_IAE": r["S8b"]["o"]["IAE"], "S9_IAE": r["S9"]["o"]["IAE"]})
    print(f"  {nom_p:16s} {J_v:6.3f}  {(','.join(non) or '-'):12s} {r['S1']['o']['depassement_pct']:9.1f} "
          f"{r['S7a']['o']['e_eff_V'] * 1e3:8.1f} {r['S8b']['o']['IAE'] * 1e3:8.2f} {r['S9']['o']['IAE'] * 1e3:8.2f}")
    L_BOB, C_CONV = L_NOM, C_NOM              # circuit nominal retabli
J_LC = [x["J"] for x in SENSIBILITE[:5]]
print(f"\n  J entre {min(J_LC):.3f} et {max(J_LC):.3f} selon L ou C a +-0.1 %.")


# %% ETAPE 5 : comparaison

titre("ETAPE 5 : Fuzzy-PID face aux references")
TOUS = {"Ziegler-Nichols": REF["Ziegler-Nichols"], "meilleur PID fige": REF["meilleur PID fige"], "Fuzzy-PID": FZ}
J_TOUS = {nom: float(np.mean(termes(res) / T_ZN)) for nom, res in TOUS.items()}
print("  " + " " * 28 + "".join(f"{nom:>20s}" for nom in TOUS))
print("  " + f"{'cout J':28s}" + "".join(f"{J_TOUS[n]:20.3f}" for n in TOUS))
print("  " + f"{'evenements non revenus':28s}" + "".join(
    f"{(','.join(c for c in CODES_IAE if not TOUS[n][c]['revenus']) or '-'):>20s}" for n in TOUS))
for code in ESSAIS:
    for cle, etiquette, fmt, k in (("depassement_pct", "depassement %", "{:.1f}", 1.0), ("IAE", "IAE mV.s", "{:.2f}", 1e3),
                                   ("e_eff_V", "ecart efficace mV", "{:.1f}", 1e3), ("iL_max_dem_A", "iL crete dem. A", "{:.1f}", 1.0)):
        if cle == "depassement_pct" and code not in ("S1", "S8a", "S8b"):
            continue
        if cle == "iL_max_dem_A" and code not in ("S1", "S8b"):
            continue
        print("  " + f"{code + ' ' + etiquette:28s}" + "".join(f"{fmt.format(TOUS[n][code]['o'][cle] * k):>20s}" for n in TOUS))
print("  " + f"{'activite |dd| (S1)':28s}" + "".join(f"{TOUS[n]['S1']['o']['dd_moyen']:20.4f}" for n in TOUS))


def decomposition(res):
    r = termes(res) / T_ZN
    return float(np.mean(r[:len(CODES_IAE)])), float(np.mean(r[len(CODES_IAE):]))


DECOMP = {nom: decomposition(TOUS[nom]) for nom in TOUS}
print("\n  Decomposition de J (moyenne des 10 IAE apres 30 ms ; moyenne des 3 IAE de demarrage) :")
for nom in TOUS:
    print(f"    {nom:20s} apres 30 ms {DECOMP[nom][0]:.3f} ; demarrage {DECOMP[nom][1]:.3f}")


# %% ETAPE 5b : essais complementaires, hors du cout J
# S10 est simule ici avec le meme regulateur, pour figurer comme les onze essais dans
# predictions_banc_fuzzy_pid.json (lu par Simuler_Fuzzy_PID.m) et dans
# banc_fuzzy_pid_resultats.json. Il n'entre ni dans J, ni dans la sensibilite, calcules plus
# haut sur les onze essais seulement.

ESSAIS_HORS_J = {sc["code"]: sc for sc in SCENARIOS_COMPLEMENTAIRES}
HORS_J = {}
if ESSAIS_HORS_J:
    titre("ETAPE 5b : essais complementaires, hors du cout J (" + ", ".join(ESSAIS_HORS_J) + ")")
    HORS_J = {"Ziegler-Nichols": evaluer(lambda: PIDClassique(), essais=ESSAIS_HORS_J),
              "Fuzzy-PID": evaluer(lambda: FuzzyPID(), essais=ESSAIS_HORS_J)}
    for code in ESSAIS_HORS_J:
        for n, res in HORS_J.items():
            print(f"  {code} {n:16s} {resume(res[code]['o'])}")
            for ev in res[code]["o"]["evenements"]:
                etat = ("reste dans la bande" if ev["reste_dans_bande"] else
                        (f"retour en {ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU"))
                print(f"      evenement a {ev['t_ms']:6.1f} ms : IAE {ev['IAE'] * 1e3:7.2f} mV.s, ecart max "
                      f"{ev['e_max_V']:6.2f} V, {etat}")


# %% ETAPE 6 : fichiers pour Simulink et MATLAB, resultats et figures

titre("ETAPE 6 : fichiers pour Simulink et MATLAB")
pas_1ms = int(round(1e-3 / TC))               # 220 periodes Tc
# 1. Resultats attendus (format des autres fichiers de predictions)
REGLAGES = {"KU": KU, "TU": TU, "KP_MIN": KP_MIN, "KP_MAX": KP_MAX, "KD_MIN": KD_MIN, "KD_MAX": KD_MAX,
            "E_M": E_M, "DE_M": DE_M, "TABLE_KP": TABLE_KP.tolist(), "TABLE_KD": TABLE_KD.tolist(),
            "TABLE_ALPHA": TABLE_ALPHA.tolist()}
predictions = {"format": 2, "Te": TC,
               "regulateur": "Fuzzy-PID (ordonnancement flou de Zhao, Tomizuka et Isaka 1993, bloc PID de "
                             "Ziegler-Nichols a gains externes)",
               "reglages": REGLAGES, "controle": {"ecart_gains_forces_PID_classique": ecart_zn}, "essais": {}}
for code, r in list(FZ.items()) + list(HORS_J.get("Fuzzy-PID", {}).items()):   # onze essais, puis S10
    K = np.array(r["reg"].K_pas)              # gains [P I D] utilises a chaque pas
    predictions["essais"][code] = {"grandeurs": r["o"],
                                   "v_toutes_les_ms": r["sim"]["v"][::pas_1ms].tolist(),
                                   "iL_toutes_les_ms": r["sim"]["iL"][::pas_1ms].tolist(),
                                   "d_moyen_par_ms": [float(np.mean(r["sim"]["d"][i:i + pas_1ms]))
                                                      for i in range(0, len(r["sim"]["d"]) - pas_1ms + 1, pas_1ms)],
                                   "K_toutes_les_ms": K[::pas_1ms].tolist()}
with open(os.path.join(DOSSIER_FUZZY, "predictions_banc_fuzzy_pid.json"), "w", encoding="utf-8") as f:
    json.dump(predictions, f, indent=1, allow_nan=False)
print("  predictions_banc_fuzzy_pid.json ecrit.")

# 2. Reglages du bloc MATLAB
savemat(os.path.join(DOSSIER_FUZZY, "fuzzy_pid_reglages.mat"),
        {"KP_MIN": KP_MIN, "KP_MAX": KP_MAX, "KD_MIN": KD_MIN, "KD_MAX": KD_MAX, "E_M": E_M, "DE_M": DE_M,
         "TABLE_KP": TABLE_KP.astype(float), "TABLE_KD": TABLE_KD.astype(float), "TABLE_ALPHA": TABLE_ALPHA,
         "K_ZN": K_ZN.reshape(1, 3), "TC": TC, "FORMAT": 2.0})
print("  fuzzy_pid_reglages.mat ecrit (reglages du bloc).")

# 3. Reference du test de rejeu MATLAB
rejeu = {}
for code in CODES_REJEU:
    r = FZ[code]
    rejeu[code] = {"e": (r["sim"]["consigne"] - r["sim"]["mesure"]).reshape(-1, 1),
                   "K": np.array(r["reg"].K_pas), "ZH": np.array(r["reg"].zh_pas)}
grille = np.array([-150, -100, -70, -50, -33.3, -10, -1, -0.05, 0, 0.05, 1, 10, 33.3, 50, 70, 100, 150], float)
grille_d = np.array([-1.0, -0.5, -0.3, -0.2, -0.1, -0.02, -0.001, 0, 0.001, 0.02, 0.1, 0.2, 0.3, 0.5, 1.0], float)
pts = np.array([[a, b] for a in grille for b in grille_d])
rejeu["FIS_ENTREES"] = pts
rejeu["FIS_SORTIES"] = np.array([zhao(a, b) for a, b in pts])
rejeu["FIS_GAINS"] = np.array([gains_zhao(a, b) for a, b in pts])
savemat(os.path.join(DOSSIER_FUZZY, "reference_rejeu_fuzzy_pid.mat"), rejeu, do_compression=True)
print(f"  reference_rejeu_fuzzy_pid.mat ecrit ({', '.join(CODES_REJEU)}, et {len(pts)} points du systeme flou).")


# 4. Resultats pour le memoire
def resume_essais(res):
    return {c: {"grandeurs": res[c]["o"], "IAE_dem": res[c]["IAE_dem"], "revenus": res[c]["revenus"]} for c in res}


sortie = {"reglages": REGLAGES, "J": J_TOUS, "decomposition": DECOMP, "sensibilite": SENSIBILITE,
          "references": {n: resume_essais(TOUS[n]) for n in TOUS},
          "essais_hors_J": {"codes": list(ESSAIS_HORS_J),
                            "note": "essais complementaires (S10), hors du cout J",
                            "references": {n: resume_essais(HORS_J[n]) for n in HORS_J}},
          "duree_s": round(time.time() - T_DEBUT, 1)}
with open(os.path.join(DOSSIER_FUZZY, "banc_fuzzy_pid_resultats.json"), "w", encoding="utf-8") as f:
    json.dump(sortie, f, indent=1, allow_nan=False, default=lambda x: None)
print("  banc_fuzzy_pid_resultats.json ecrit.")

# 5. Figure : tension, gains et courant sur S1 (demarrage), S8b et S9
fig, axes = plt.subplots(3, 3, figsize=(16, 9), sharex="col")
for col, (code, xlim) in enumerate((("S1", (0, 15)), ("S8b", None), ("S9", None))):
    for nom, couleur in (("Ziegler-Nichols", "0.55"), ("meilleur PID fige", "tab:green"), ("Fuzzy-PID", "tab:purple")):
        sim = TOUS[nom][code]["sim"]
        axes[0, col].plot(sim["t"] * 1e3, sim["v"], lw=0.6, color=couleur, label=nom)
        axes[2, col].plot(sim["t"] * 1e3, sim["iL"], lw=0.4, color=couleur)
    axes[0, col].plot(sim["t"] * 1e3, sim["consigne"], "k--", lw=0.6)
    K = np.array(FZ[code]["reg"].K_pas) / K_ZN
    t_ms = FZ[code]["sim"]["t"] * 1e3
    for i, nom in enumerate(("P", "I", "D")):
        axes[1, col].plot(t_ms, K[:, i], lw=0.6, label=nom)
    axes[0, col].set_ylim(30, 200)
    axes[0, col].set_title(f"{code} : tension de sortie", fontsize=9, loc="left")
    axes[1, col].set_title(f"{code} : gains du Fuzzy-PID (multiples de Ziegler-Nichols)", fontsize=9, loc="left")
    axes[2, col].set_title(f"{code} : courant de la bobine (observe)", fontsize=9, loc="left")
    axes[2, col].set_xlabel("temps (ms)")
    if xlim:
        axes[2, col].set_xlim(*xlim)
    for ax in axes[:, col]:
        ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7)
axes[1, 0].legend(fontsize=7)
axes[0, 0].set_ylabel("V")
axes[2, 0].set_ylabel("A")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_FUZZY, "banc_fuzzy_pid.png"), dpi=110)
plt.close(fig)
print("  banc_fuzzy_pid.png ecrit.")

# 6. Figure : K'p, K'd et alpha en fonction de e' et De' (pour le memoire)
xv = np.linspace(-1.2, 1.2, 241)
Z = np.array([[zhao(a * E_M, b * DE_M) for b in xv] for a in xv])   # 241 x 241 x 3
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for i, (nom, ax) in enumerate(zip(("K'p", "K'd", "alpha"), axes)):
    im = ax.pcolormesh(xv, xv, Z[:, :, i], cmap="viridis", shading="nearest")
    ax.set_xlabel("De' = De / DE_M")
    ax.set_ylabel("e' = e / E_M")
    ax.set_title(nom, fontsize=9, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_FUZZY, "fuzzy_pid_surfaces.png"), dpi=110)
plt.close(fig)
print(f"  fuzzy_pid_surfaces.png ecrit. Duree totale : {time.time() - T_DEBUT:.0f} s.")
