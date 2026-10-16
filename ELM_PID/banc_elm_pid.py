# =============================================================================
# banc_elm_pid.py
#
# OBJECTIF
#   Faire tourner l'ELM-PID sur le banc commun, sur les onze essais de
#   developpement et sur S10, a cote des references (PID de Ziegler-Nichols,
#   qui est aussi l'ELM-PID sans adaptation, et meilleur PID fige sur
#   grille), mesurer ce que chaque piece apporte (ablations) et ecrire ce
#   que Simulink doit retrouver :
#     - predictions_banc_elm_pid.json (lu par Simuler_ELM_PID.m) ;
#     - elm_pid_reglages.mat (reglages lus par elm_pid_adaptatif.m) ;
#     - reference_rejeu_elm_pid.mat (lu par Tester_ELM_PID_Rejeu.m).
#
# L'ELM-PID, TEL QU'IL EST CODE ICI ET DANS elm_pid_adaptatif.m
#   A chaque periode Tc = 1/220 000 s :
#   1. Le bloc PID de la base commune (Parallel, derivee filtree N =
#      64 122.9 rad/s, Forward Euler, sortie dans [0.01 ; 0.99], clamping,
#      integrateur a 0.5 et filtre a 0.01 au depart) calcule u a partir de
#      e = consigne - mesure, avec les gains P, I, D que lui donne
#      l'adaptateur. Sans adaptation, c'est le PID de Ziegler-Nichols.
#   2. L'adaptateur recoit e, la mesure de Vout (jamais 100 - e) et u. Sur
#      des fenetres de 0.5 ms (110 x Tc) : moyennes ybar, dbar, ebar, et
#      sensibilite moyenne s de u aux trois gains (derivee exacte de la loi
#      du bloc depuis le debut de la fenetre). Une fenetre ou u est en
#      butee est suspecte.
#   3. En fin de fenetre : l'ELM (elm_pid_modele.mat) predit ybar(n) a
#      partir de [ybar(n-1), ybar(n-2), dbar(n), dbar(n-1), dbar(n-2)] et
#      donne J = d ybar(n) / d dbar(n).
#   4. Porte : ouverte si ni la fenetre ni les deux precedentes ne sont
#      suspectes. Porte ouverte :
#      a. OS-ELM (moindres carres recursifs, sans oubli) si
#         |ybar - prediction| > ZONE_MORTE_ESTIMATION ;
#      b. si |ebar| > ZONE_MORTE et J > 0 : multiplicateurs de
#         Ziegler-Nichols x = K / K_ZN,
#           phi = J (s .* K_ZN),
#           dx = eta ebar phi / (eps + phi.phi) + alpha (x - x_prec),
#         puis projection sur l'ensemble admissible
#         (ensemble_gains_elm.mat).
#   Les gains sortis par l'adaptateur ne dependent que de son etat : ceux
#   d'une fin de fenetre servent a partir du pas suivant.
#
# COMMENT LANCER CE SCRIPT
#   python banc_elm_pid.py      (une a deux minutes : ablations et
#   sensibilite comprises)
#
# ORDRE D'EXECUTION
#   1. entrainement_elm.py ; 2. ensemble_gains_elm.py ; 3. ce script ;
#   4. Tester_ELM_PID_Rejeu.m ; 5. Construction_ELM_PID.m ;
#   6. Simuler_ELM_PID.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques, banc commun, modele, ensemble admissible

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
with open(os.path.join(DOSSIER_ELM, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC or "iL_max_dem_A" not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py ne contient pas les definitions attendues (marque de fin, observation du courant).")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])

T_DEBUT = time.time()
ESSAIS = {sc["code"]: sc for sc in SCENARIOS}
CODES_IAE = ["S1", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8a", "S8b", "S9"]
CODES_DEM = ["S1", "S8a", "S8b"]
K30 = int(round(0.03 / TC))
NF = 110
CODES_REJEU = ["S1", "S3", "S7b", "S8b", "S9"]
K_ZN = np.array([PID_P, PID_I, PID_D])   # Ziegler-Nichols, gains du bloc PID [P, I, D]

REGLAGES = {"ETA": 0.5, "ZONE_MORTE": 0.1, "ZONE_MORTE_ESTIMATION": 0.1, "ALPHA": 0.001, "EPS_PHI": 1e-3,
            "LAMBDA": 1.0, "APPRENDRE": 1.0, "GLISSEMENT": 1.0, "JACOBIEN_CONSTANT": float("nan"),
            "MULT_MIN": 0.25, "MULT_MAX": 4.0, "N_BISSECTIONS": 30.0, "H_GRADIENT": 1e-4,
            "RESTAURATION": 1.0, "N_PROJECTION": 10.0}

with open(os.path.join(DOSSIER_ELM, "predictions_banc_pid_fige.json"), encoding="utf-8") as f:
    GAINS_FIGE = json.load(f)["gains"]
with open(os.path.join(DOSSIER_ELM, "entrainement_elm_resultats.json"), encoding="utf-8") as f:
    APPR = json.load(f)
J_REGRESSION = float(APPR["references"]["regression lineaire, ordre 2"]["J_median"])   # constante

MODELE = loadmat(os.path.join(DOSSIER_ELM, "elm_pid_modele.mat"))
if int(MODELE["NF"].item()) != NF or int(MODELE["ORDRE"].item()) != 2:
    raise RuntimeError("elm_pid_modele.mat n'est pas le modele attendu (fenetres de 110 x Tc, ordre 2).")
ENS = loadmat(os.path.join(DOSSIER_ELM, "ensemble_gains_elm.mat"))
G_TABLE = ENS["G_TABLE"].astype(float)
LOG2_MIN, PAS_LOG2 = float(ENS["LOG2_MIN"].item()), float(ENS["PAS_LOG2"].item())
N_TABLE = G_TABLE.shape[0]
if np.max(np.abs(ENS["K_DEPART"].ravel() - K_ZN) / K_ZN) > 1e-12:
    raise RuntimeError("ensemble_gains_elm.mat n'a pas ete calcule avec les gains de Ziegler-Nichols de ce banc.")
print(f"ELM : {MODELE['W'].shape[1]} neurones, ordre {int(MODELE['ORDRE'].item())} ; ensemble admissible : "
      f"{N_TABLE}^3 points, {np.mean(G_TABLE <= 0):.1%} admissibles ; jacobien constant des ablations "
      f"{J_REGRESSION:.4f}.")


def titre(texte):
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78)


# %% ETAPE 1 : bloc PID, projection et adaptateur ELM (copie de elm_pid_adaptatif.m)

TCN = TC * PID_N


class BlocPID:
    """Le bloc "PID Controller" de la base commune (memes operations que
    PIDClassique du banc commun), gains P, I, D recus de l'exterieur a chaque
    pas (ports de gains externes du bloc Simulink, comme pour le Fuzzy-PID).
    Conditions initiales du bloc : integrateur 0.5, filtre 0.01."""

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


class AdaptateurELM:
    """Copie Python de elm_pid_adaptatif.m (memes reglages, memes noms,
    meme ordre des operations : Tester_ELM_PID_Rejeu.m compare les deux pas
    par pas). Entrees a chaque pas : erreur e, mesure de Vout, commande u
    sortie du bloc PID. Sortie : gains [P I D] du bloc PID. Les gains sortis
    au pas k ne dependent que de l'etat (aucune traversee directe) : ils
    changent en fin de fenetre et servent a partir du pas suivant. Les
    arguments servent aux ablations."""

    def __init__(self, **modif):
        self.r = dict(REGLAGES, **modif)
        r = self.r
        self.eta, self.zm, self.zme = r["ETA"], r["ZONE_MORTE"], r["ZONE_MORTE_ESTIMATION"]
        self.alpha, self.eps, self.lam = r["ALPHA"], r["EPS_PHI"], r["LAMBDA"]
        self.apprendre, self.Jc = bool(r["APPRENDRE"]), r["JACOBIEN_CONSTANT"]
        self.adapter = r.get("ADAPTER", True)
        self.W = MODELE["W"].astype(float)
        self.B = MODELE["B"].astype(float).ravel()
        self.BETA = MODELE["BETA"].astype(float).ravel().copy()
        self.P = MODELE["P0"].astype(float).copy()
        self.XM = MODELE["X_MOY"].astype(float).ravel()
        self.XE = MODELE["X_EC"].astype(float).ravel()
        self.TM = float(MODELE["T_MOY"].item())
        self.TE_ = float(MODELE["T_EC"].item())
        self.n = self.W.shape[1]
        self.x = np.ones(3)                   # multiplicateurs de Ziegler-Nichols
        self.x_prec = self.x.copy()
        self.K = K_ZN * self.x
        self.cpt, self.sy, self.sd, self.se, self.hors = 0, 0.0, 0.0, 0.0, False
        self.SE, self.phiD = 0.0, 0.0          # somme des erreurs et filtre de sensibilite, remis a zero par fenetre
        self.sS1 = self.sS2 = self.sS3 = 0.0
        self.yb1 = self.yb2 = self.db1 = self.db2 = 0.0
        self.sus1 = self.sus2 = True
        self.nfen = 0
        self.journal = []

    def modele(self, x):
        xn = (x - self.XM) / self.XE
        gg = 1.0 / (1.0 + np.exp(-(xn @ self.W + self.B)))
        h = np.concatenate([gg, xn, [1.0]])
        y = (h @ self.BETA) * self.TE_ + self.TM
        dy = ((gg * (1.0 - gg)) * self.BETA[:self.n]) @ self.W.T + self.BETA[self.n:self.n + 5]
        return y, dy[2] / self.XE[2] * self.TE_, h

    def mise_a_jour(self, e, mesure, u):
        # Sensibilite de u(k) a un changement des gains au debut de la fenetre
        # (derivee exacte de la loi du bloc, hors clamping) :
        #   du/dP = e(k) ; du/dI = Tc x (somme des e de la fenetre avant k) ;
        #   du/dD = N (e(k) - phiD(k)), phiD(k+1) = phiD(k) + Tc N (e(k) - phiD(k)), phiD = 0 au debut.
        self.sS1 += e
        self.sS2 += TC * self.SE
        self.sS3 += PID_N * (e - self.phiD)
        self.SE += e
        self.phiD += TCN * (e - self.phiD)
        self.sy, self.sd, self.se = self.sy + mesure, self.sd + u, self.se + e
        if u <= D_MIN or u >= D_MAX:          # sortie du bloc en butee : la fenetre est suspecte
            self.hors = True
        self.cpt += 1
        if self.cpt >= NF:
            self.fin_de_fenetre()

    def fin_de_fenetre(self):
        yb, db, eb = self.sy / NF, self.sd / NF, self.se / NF
        s = np.array([self.sS1, self.sS2, self.sS3]) / NF
        sus, porte, err, J, adapte, appris, change = self.hors, False, float("nan"), float("nan"), False, False, False
        if self.nfen >= 2:
            xin = np.array([self.yb1, self.yb2, db, self.db1, self.db2])
            yhat, J, h = self.modele(xin)
            if not np.isnan(self.Jc):
                J = self.Jc
            err = yb - yhat
            porte = (not sus) and (not self.sus1) and (not self.sus2)
            if porte:
                if self.apprendre and abs(err) > self.zme:          # a. OS-ELM, zone morte d'estimation
                    Ph = self.P @ h
                    gain = Ph / (self.lam + h @ Ph)
                    self.BETA = self.BETA + gain * (err / self.TE_)
                    self.P = (self.P - np.outer(gain, Ph)) / self.lam
                    appris = True
                ancien = self.x.copy()
                if self.adapter and abs(eb) > self.zm and J > 0:     # b. gains
                    phi = J * (s * K_ZN)
                    dx = self.eta * eb * phi / (self.eps + phi @ phi) + self.alpha * (self.x - self.x_prec)
                    self.x = projeter(self.x, dx, self.r)
                    adapte = True
                    change = bool(np.any(self.x != ancien))
                self.x_prec = ancien
                self.K = K_ZN * self.x
        self.journal.append({"ybar": yb, "dbar": db, "ebar": eb, "err": err, "J": J, "porte": porte,
                             "adapte": adapte, "appris": appris, "change": change, "x": self.x.copy(),
                             "K": self.K.copy(), "s": s})
        self.yb2, self.yb1 = self.yb1, yb
        self.db2, self.db1 = self.db1, db
        self.sus2, self.sus1 = self.sus1, sus
        self.nfen += 1
        self.cpt, self.sy, self.sd, self.se, self.hors = 0, 0.0, 0.0, 0.0, False
        self.SE, self.phiD = 0.0, 0.0
        self.sS1 = self.sS2 = self.sS3 = 0.0


class ELMPID:
    """Le montage Simulink : adaptateur ELM -> gains -> bloc
    PID de la base commune -> u, et u renvoye a l'adaptateur."""
    utilise_mesure = True

    def __init__(self, **modif):
        self.adapt = AdaptateurELM(**modif)
        self.bloc = BlocPID()
        self.journal = self.adapt.journal
        self.K_pas = []

    @property
    def K(self):
        return self.adapt.K

    def pas(self, e, mesure):
        K = self.adapt.K                       # sortie de l'adaptateur (son etat seulement)
        self.K_pas.append(K)
        u = self.bloc.pas(e, K)
        self.adapt.mise_a_jour(e, mesure, u)
        return u


essai = ELMPID()
ecart_y = max(abs(essai.adapt.modele(MODELE["X_TEST"][i].astype(float))[0] - MODELE["Y_TEST"][i, 0])
              for i in range(MODELE["X_TEST"].shape[0]))
ecart_J = max(abs(essai.adapt.modele(MODELE["X_TEST"][i].astype(float))[1] - MODELE["J_TEST"][i, 0])
              for i in range(MODELE["X_TEST"].shape[0]))
print(f"Vecteurs de test du modele : ecart {ecart_y:.1e} V sur ybar, {ecart_J:.1e} sur le jacobien.")
if ecart_y > 1e-9 or ecart_J > 1e-9:
    raise RuntimeError("elm_pid_modele.mat ne redonne pas ses vecteurs de test.")
if g_interp(np.ones(3)) > 0.0:
    raise RuntimeError("Le depart (gains de Ziegler-Nichols) n'est pas admissible dans la table.")


# %% ETAPE 2 : evaluation et references

def evaluer(fabrique, codes=None, essais=None):
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
    return np.array([res[c]["o"]["IAE"] for c in CODES_IAE] + [res[c]["IAE_dem"] for c in CODES_DEM])


class PIDParallele(PIDClassique):
    """PID du bloc Simulink (forme Parallel, clamping) avec d'autres gains."""

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


def fenetres_adaptees(reg):
    """Fenetres ou un pas de gradient a ete calcule (porte ouverte, |ebar| > zone morte, J > 0)."""
    return [i for i, x in enumerate(reg.journal) if x["adapte"]]


def fenetres_changees(reg):
    """Fenetres ou les gains ont effectivement change apres projection."""
    return [i for i, x in enumerate(reg.journal) if x["change"]]


def bilan(nom, res):
    J = float(np.mean(termes(res) / T_ZN))
    non = [c for c in CODES_IAE if not res[c]["revenus"]]
    print(f"\n  {nom} : J = {J:.3f} ; {'tous les evenements reviennent' if not non else 'pas revenus : ' + ', '.join(non)}")
    for code, r in res.items():
        o = r["o"]
        ligne = (f"    {code:4s} depassement {o['depassement_pct']:5.1f} %, IAE {o['IAE'] * 1e3:7.2f} mV.s, "
                 f"efficace {o['e_eff_V'] * 1e3:7.1f} mV, |dd| {o['dd_moyen']:.4f}")
        if hasattr(r["reg"], "journal"):
            jr = r["reg"].journal
            Kf = jr[-1]["K"] / K_ZN
            ad = fenetres_adaptees(r["reg"])
            ch = fenetres_changees(r["reg"])
            ligne += (f" ; porte {100 * np.mean([x['porte'] for x in jr]):3.0f} %, pas calcules {len(ad):2d}, gains "
                      f"changes {len(ch):2d} (dont {sum(1 for i in ch if (i + 1) * NF > K30):2d} apres 30 ms), finaux "
                      f"x({Kf[0]:.3f}, {Kf[1]:.3f}, {Kf[2]:.3f})")
        print(ligne)
    return J


titre("ETAPE 2 : les references")
REF = {"Ziegler-Nichols": evaluer(lambda: PIDClassique()),
       "meilleur PID fige": evaluer(lambda: PIDParallele(GAINS_FIGE["P"], GAINS_FIGE["I"], GAINS_FIGE["D"]))}
T_ZN = termes(REF["Ziegler-Nichols"])
J_FIGE = float(np.mean(termes(REF["meilleur PID fige"]) / T_ZN))
print(f"  Meilleur PID fige sur grille : J = {J_FIGE:.3f}.")


# %% ETAPE 3 : ELM-PID sur les onze essais

titre("ETAPE 3 : ELM-PID (le point de depart, gains figes, est le PID de Ziegler-Nichols lui-meme)")
ELM = evaluer(lambda: ELMPID())
J_ELM = bilan("ELM-PID", ELM)
fige = evaluer(lambda: ELMPID(ADAPTER=False))
ecart_zn = max(float(np.max(np.abs(fige[c]["sim"]["d"] - REF["Ziegler-Nichols"][c]["sim"]["d"]))) for c in fige)
print(f"\n  Controle : ELM-PID sans adaptation contre le PID de Ziegler-Nichols, onze essais, ecart maximal sur d : {ecart_zn:.1e}.")
if ecart_zn != 0.0:
    raise RuntimeError("Sans adaptation, l'ELM-PID ne redonne pas le PID de Ziegler-Nichols a l'identique.")
for code in ("S1", "S2", "S3"):
    reg = ELM[code]["reg"]
    for i in fenetres_adaptees(reg):
        x = reg.journal[i]
        print(f"  {code} : pas calcule a t = {(i + 1) * 0.5:5.1f} ms, ebar {x['ebar']:+.3f} V, jacobien du reseau "
              f"{x['J']:.2f} ; gains {'changes' if x['change'] else 'inchanges'} -> x({x['x'][0]:.3f}, {x['x'][1]:.3f}, "
              f"{x['x'][2]:.3f})")
admis = [g_interp(x["x"]) <= 0.0 for c in ELM for x in ELM[c]["reg"].journal]
print(f"  Controle : gains admissibles a la fin de toutes les fenetres des onze essais : {'oui' if all(admis) else 'NON'}.")


# %% ETAPE 4 : ablations

titre("ETAPE 4 : ce que chaque piece apporte (ablations)")
ABLATIONS = {"Jacobien constant": {"JACOBIEN_CONSTANT": J_REGRESSION},
             "Sans OS-ELM": {"APPRENDRE": 0.0},
             "Sans restauration": {"RESTAURATION": 0.0},
             "Sans glissement": {"GLISSEMENT": 0.0}}
RES_ABL, J_ABL = {}, {}
for nom, modif in ABLATIONS.items():
    RES_ABL[nom] = evaluer(lambda: ELMPID(**modif))
    J_ABL[nom] = float(np.mean(termes(RES_ABL[nom]) / T_ZN))
print(f"\n  {'variante':28s} {'J':>6s}  {'non revenus':12s} {'S3':>7s} {'S5':>7s} {'S8b':>7s} {'S9':>8s} {'gains changes':>13s}")
print("  (IAE de 30 ms a la fin en mV.s ; gains changes : fenetres ou les gains ont change, total des onze essais)")
for nom, res, J in ([("ELM-PID", ELM, J_ELM)] + [(n, RES_ABL[n], J_ABL[n]) for n in ABLATIONS]):
    non = [c for c in CODES_IAE if not res[c]["revenus"]]
    nb = sum(len(fenetres_changees(res[c]["reg"])) for c in res) if hasattr(res["S1"]["reg"], "journal") else 0
    print(f"  {nom:28s} {J:6.3f}  {(','.join(non) or '-'):12s} {res['S3']['o']['IAE'] * 1e3:7.2f} "
          f"{res['S5']['o']['IAE'] * 1e3:7.2f} {res['S8b']['o']['IAE'] * 1e3:7.2f} {res['S9']['o']['IAE'] * 1e3:8.2f} {nb:13d}")


# %% ETAPE 4b : sensibilite a de tres petits ecarts du circuit

titre("ETAPE 4b : sensibilite (L ou C a +-0.1 %)")
L_NOM, C_NOM = L_BOB, C_CONV
SENSIBILITE = []
for nom_p, fl, fc in (("nominal", 1.0, 1.0), ("L - 0.1 %", 0.999, 1.0), ("L + 0.1 %", 1.001, 1.0),
                      ("C - 0.1 %", 1.0, 0.999), ("C + 0.1 %", 1.0, 1.001)):
    L_BOB, C_CONV = L_NOM * fl, C_NOM * fc
    r_e = ELM if nom_p == "nominal" else evaluer(lambda: ELMPID())
    J_e = float(np.mean(termes(r_e) / T_ZN))
    rev = all(r_e[c]["revenus"] for c in CODES_IAE)
    nb = sum(len(fenetres_changees(r_e[c]["reg"])) for c in r_e)
    SENSIBILITE.append({"circuit": nom_p, "J": J_e, "tous_revenus": rev, "fenetres_gains_changes": nb})
    print(f"  {nom_p:12s} J = {J_e:.3f} ; gains changes sur {nb} fenetres ; {'tous revenus' if rev else 'PAS TOUS REVENUS'}")
L_BOB, C_CONV = L_NOM, C_NOM


# %% ETAPE 5 : comparaison

titre("ETAPE 5 : ELM-PID face aux references")
TOUS = {"Ziegler-Nichols": REF["Ziegler-Nichols"], "meilleur PID fige": REF["meilleur PID fige"], "ELM-PID": ELM}
J_TOUS = {n: float(np.mean(termes(r) / T_ZN)) for n, r in TOUS.items()}
print("  " + " " * 28 + "".join(f"{n:>20s}" for n in TOUS))
print("  " + f"{'cout J':28s}" + "".join(f"{J_TOUS[n]:20.3f}" for n in TOUS))
print("  " + f"{'evenements non revenus':28s}" + "".join(
    f"{(','.join(c for c in CODES_IAE if not TOUS[n][c]['revenus']) or '-'):>20s}" for n in TOUS))
for code in ESSAIS:
    for cle, etiquette, fmt, k in (("depassement_pct", "depassement %", "{:.1f}", 1.0), ("IAE", "IAE mV.s", "{:.2f}", 1e3),
                                   ("e_max_V", "ecart maximal V", "{:.2f}", 1.0), ("e_eff_V", "ecart efficace mV", "{:.1f}", 1e3),
                                   ("iL_max_dem_A", "iL crete dem. A", "{:.1f}", 1.0)):
        if cle == "depassement_pct" and code not in ("S1", "S8a", "S8b"):
            continue
        if cle == "iL_max_dem_A" and code not in ("S1", "S8b"):
            continue
        if cle == "e_max_V" and code not in ("S2", "S3"):
            continue
        print("  " + f"{code + ' ' + etiquette:28s}" + "".join(f"{fmt.format(TOUS[n][code]['o'][cle] * k):>20s}" for n in TOUS))
DECOMP = {}
for n in TOUS:
    r = termes(TOUS[n]) / T_ZN
    DECOMP[n] = (float(np.mean(r[:10])), float(np.mean(r[10:])))
print("  Decomposition de J (apres 30 ms ; demarrage) : " +
      " ; ".join(f"{n} {DECOMP[n][0]:.3f} / {DECOMP[n][1]:.3f}" for n in TOUS))



# %% ETAPE 5b : essais complementaires, hors du cout J
# S10 est simule ici avec le meme regulateur, pour figurer comme les onze essais dans
# predictions_banc_elm_pid.json (lu par Simuler_ELM_PID.m) et dans banc_elm_pid_resultats.json.
# Il n'entre ni dans J, ni dans les ablations ou la sensibilite, calcules plus haut sur ESSAIS
# (les onze essais) seulement.

ESSAIS_HORS_J = {sc["code"]: sc for sc in SCENARIOS_COMPLEMENTAIRES}
HORS_J = {}
if ESSAIS_HORS_J:
    titre("ETAPE 5b : essais complementaires, hors du cout J (" + ", ".join(ESSAIS_HORS_J) + ")")
    HORS_J = {"Ziegler-Nichols": evaluer(lambda: PIDClassique(), essais=ESSAIS_HORS_J),
              "ELM-PID": evaluer(lambda: ELMPID(), essais=ESSAIS_HORS_J)}
    for code in ESSAIS_HORS_J:
        for n, res in HORS_J.items():
            print(f"  {code} {n:16s} {resume(res[code]['o'])}")
            for ev in res[code]["o"]["evenements"]:
                etat = ("reste dans la bande" if ev["reste_dans_bande"] else
                        (f"retour en {ev['t_retour_ms']:.2f} ms" if ev["revenu"] else "PAS REVENU"))
                print(f"      evenement a {ev['t_ms']:6.1f} ms : IAE {ev['IAE'] * 1e3:7.2f} mV.s, ecart max "
                      f"{ev['e_max_V']:6.2f} V, {etat}")
        reg = HORS_J["ELM-PID"][code]["reg"]
        Kf = reg.journal[-1]["K"] / K_ZN
        print(f"  {code} ELM-PID : gains changes sur {len(fenetres_changees(reg))} fenetres ; finaux "
              f"x({Kf[0]:.3f}, {Kf[1]:.3f}, {Kf[2]:.3f})")


# %% ETAPE 6 : fichiers pour Simulink et MATLAB, resultats et figure

titre("ETAPE 6 : fichiers pour Simulink et MATLAB")
pas_1ms = int(round(1e-3 / TC))
REGL_JSON = {k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in REGLAGES.items()}
predictions = {"format": 2, "Te": TC,               # format des predictions (lu par Simuler_ELM_PID.m)
               "regulateur": "ELM-PID (bloc PID de la base commune a gains externes, ELM de Lu adapte, porte, projection)",
               "reglages": REGL_JSON, "K_depart": K_ZN.tolist(),
               "controle": {"ecart_vecteurs_test_V": ecart_y, "ecart_sans_adaptation_ZN": ecart_zn},
               "essais": {}}
for code, r in list(ELM.items()) + list(HORS_J.get("ELM-PID", {}).items()):   # onze essais, puis S10
    Kp = np.array(r["reg"].K_pas)
    predictions["essais"][code] = {"grandeurs": r["o"],
                                   "v_toutes_les_ms": r["sim"]["v"][::pas_1ms].tolist(),
                                   "iL_toutes_les_ms": r["sim"]["iL"][::pas_1ms].tolist(),
                                   "d_moyen_par_ms": [float(np.mean(r["sim"]["d"][i:i + pas_1ms]))
                                                      for i in range(0, len(r["sim"]["d"]) - pas_1ms + 1, pas_1ms)],
                                   "K_toutes_les_ms": Kp[::pas_1ms].tolist(),
                                   "fenetres_adaptees": len(fenetres_adaptees(r["reg"])),
                                   "fenetres_gains_changes": len(fenetres_changees(r["reg"])),
                                   "K_final": r["reg"].journal[-1]["K"].tolist()}
with open(os.path.join(DOSSIER_ELM, "predictions_banc_elm_pid.json"), "w", encoding="utf-8") as f:
    json.dump(predictions, f, indent=1, allow_nan=False)
print("  predictions_banc_elm_pid.json ecrit.")

savemat(os.path.join(DOSSIER_ELM, "elm_pid_reglages.mat"),
        {"K_DEPART": K_ZN.reshape(1, 3), "ETA": REGLAGES["ETA"], "ZONE_MORTE": REGLAGES["ZONE_MORTE"],
         "ZONE_MORTE_ESTIMATION": REGLAGES["ZONE_MORTE_ESTIMATION"], "ALPHA": REGLAGES["ALPHA"],
         "EPS_PHI": REGLAGES["EPS_PHI"], "LAMBDA": REGLAGES["LAMBDA"], "APPRENDRE": REGLAGES["APPRENDRE"],
         "GLISSEMENT": REGLAGES["GLISSEMENT"], "JACOBIEN_CONSTANT": REGLAGES["JACOBIEN_CONSTANT"],
         "MULT_MIN": REGLAGES["MULT_MIN"], "MULT_MAX": REGLAGES["MULT_MAX"],
         "N_BISSECTIONS": REGLAGES["N_BISSECTIONS"], "H_GRADIENT": REGLAGES["H_GRADIENT"],
         "RESTAURATION": REGLAGES["RESTAURATION"], "N_PROJECTION": REGLAGES["N_PROJECTION"],
         "NF": float(NF), "TC": TC, "N_FILTRE": PID_N, "D_MIN": D_MIN, "D_MAX": D_MAX, "FORMAT": 4.0})
print("  elm_pid_reglages.mat ecrit.")

rejeu = {}
for code in CODES_REJEU:
    r = ELM[code]
    jr = r["reg"].journal
    rejeu[code] = {"e": (r["sim"]["consigne"] - r["sim"]["mesure"]).reshape(-1, 1),
                   "mesure": r["sim"]["mesure"].reshape(-1, 1), "u": r["sim"]["d"].reshape(-1, 1),
                   "K": np.array(r["reg"].K_pas), "porte": np.array([x["porte"] for x in jr], float).reshape(-1, 1),
                   "adapte": np.array([x["adapte"] for x in jr], float).reshape(-1, 1)}
savemat(os.path.join(DOSSIER_ELM, "reference_rejeu_elm_pid.mat"), rejeu, do_compression=True)
print(f"  reference_rejeu_elm_pid.mat ecrit ({', '.join(CODES_REJEU)}).")


def resume_essais(res):
    return {c: {"grandeurs": res[c]["o"], "IAE_dem": res[c]["IAE_dem"], "revenus": res[c]["revenus"]} for c in res}


sortie = {"reglages": REGL_JSON, "J": J_TOUS, "decomposition": DECOMP, "J_ablations": J_ABL,
          "J_regression": J_REGRESSION, "sensibilite": SENSIBILITE,
          "references": {n: resume_essais(TOUS[n]) for n in TOUS},
          "ablations": {n: resume_essais(RES_ABL[n]) for n in RES_ABL},
          "gains_elm": {c: [x["K"].tolist() for x in ELM[c]["reg"].journal] for c in ELM},
          "fenetres_adaptees": {c: fenetres_adaptees(ELM[c]["reg"]) for c in ELM},
          "fenetres_gains_changes": {c: fenetres_changees(ELM[c]["reg"]) for c in ELM},
          "porte_elm": {c: float(np.mean([x["porte"] for x in ELM[c]["reg"].journal])) for c in ELM},
          "essais_hors_J": {"codes": list(ESSAIS_HORS_J),
                            "note": "essais complementaires (S10), hors du cout J",
                            "references": {n: resume_essais(HORS_J[n]) for n in HORS_J},
                            "gains_elm": {c: [x["K"].tolist() for x in HORS_J["ELM-PID"][c]["reg"].journal]
                                          for c in ESSAIS_HORS_J},
                            "fenetres_adaptees": {c: fenetres_adaptees(HORS_J["ELM-PID"][c]["reg"]) for c in ESSAIS_HORS_J},
                            "fenetres_gains_changes": {c: fenetres_changees(HORS_J["ELM-PID"][c]["reg"])
                                                       for c in ESSAIS_HORS_J}},
          "duree_s": round(time.time() - T_DEBUT, 1)}
with open(os.path.join(DOSSIER_ELM, "banc_elm_pid_resultats.json"), "w", encoding="utf-8") as f:
    json.dump(sortie, f, indent=1, allow_nan=False, default=lambda x: None)
print("  banc_elm_pid_resultats.json ecrit.")

fig, axes = plt.subplots(3, 4, figsize=(20, 9.5), sharex="col")
for col, (code, xlim) in enumerate((("S1", (0, 30)), ("S2", (40, 90)), ("S3", (40, 90)), ("S8b", None))):
    for nom, couleur in (("Ziegler-Nichols", "0.6"), ("meilleur PID fige", "tab:green"),
                         ("ELM-PID", "tab:blue")):
        sim = TOUS[nom][code]["sim"]
        axes[0, col].plot(sim["t"] * 1e3, sim["v"], lw=0.6, color=couleur, label=nom)
        axes[2, col].plot(sim["t"] * 1e3, sim["iL"], lw=0.4, color=couleur)
    axes[0, col].plot(sim["t"] * 1e3, sim["consigne"], "k--", lw=0.6)
    Kp = np.array(ELM[code]["reg"].K_pas) / K_ZN
    t_ms = ELM[code]["sim"]["t"] * 1e3
    for i, nom in enumerate(("Kp", "Ki", "Kd")):
        axes[1, col].plot(t_ms, Kp[:, i], lw=1.0, label=nom)
    axes[0, col].set_ylim((85, 115) if code in ("S2", "S3") else (30, 200))
    axes[0, col].set_title(f"{code} : tension de sortie", fontsize=9, loc="left")
    axes[1, col].set_title(f"{code} : gains de l'ELM-PID (multiples de ZN)", fontsize=9, loc="left")
    axes[2, col].set_title(f"{code} : courant de la bobine (observe)", fontsize=9, loc="left")
    axes[2, col].set_xlabel("temps (ms)")
    if xlim:
        axes[2, col].set_xlim(*xlim)
    for ax in axes[:, col]:
        ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7)
axes[1, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_ELM, "banc_elm_pid.png"), dpi=110)
plt.close(fig)
print(f"  banc_elm_pid.png ecrit. Duree totale : {time.time() - T_DEBUT:.0f} s.")
