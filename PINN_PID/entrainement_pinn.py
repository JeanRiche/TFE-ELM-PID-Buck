# =============================================================================
# entrainement_pinn.py
#
# OBJECTIF
#   Apprendre le modele interne du PINN-PID : un reseau de neurones informe
#   par la physique (PINN) qui predit l'etat du convertisseur Buck sur une
#   periode du regulateur Tc = 1/220 000 s. Le PINN apprend a approcher le
#   modele a partir de sa physique et des donnees ; au mieux, il l'egale. Le
#   modele physique sert donc de reference : on mesure de combien le PINN
#   s'en ecarte.
#
# LE PINN (Ito et Wasa, arXiv:2510.04591, 2025, sections II-B et III-B)
#   Etat x = [Vout, iL]. Pendant une tranche Tc, le MOSFET conduit une
#   fraction s du temps (deduite du rapport cyclique et de la phase de la
#   porteuse, memes conventions que le banc commun). Le reseau N (tanh,
#   couches cachees de meme largeur) donne l'etat a l'instant tau de la
#   tranche, avec la condition initiale imposee par construction :
#     x(tau) = x0 + (tau / Tc) . S . N(z),
#     S = [Tc x 10 / C, Tc x 100 / L] (echelles : 10 A dans le condensateur,
#     100 V sur la bobine),
#     z = entrees normalisees [tau, Vout0, iL0, s, Vin, G] (G = 1/R).
#   Perte L = L_donnees + lambda L_physique, lambda = 1 (valeur d'Ito et
#   Wasa) :
#     L_donnees : ecart, a tau = Tc, entre la prediction et l'etat enregistre
#       au pas suivant (transitions des enregistrements Simulink ELM_A1 et
#       ELM_A2), divise par S ;
#     L_physique : residu (dx/dt - f(x)) Tc / S de l'equation du modele
#       moyen par tranche, en des points tires dans le domaine :
#         C dv/dt = iL - G v,
#         L diL/dt = s Vin - v - (1 - s) Vf - (s Ron + (1 - s) Rd) iL,
#         diL/dt = 0 si iL <= 0 et diL/dt < 0 (blocage de la diode).
#   Optimisation : L-BFGS (Ito et Wasa), 15 000 iterations au plus.
#   Gradients calcules a la main (retropropagation, avec la derivee par
#   rapport a tau propagee en mode direct pour le residu physique) ; ils
#   sont controles par differences finies au debut.
#
# DOMAINE D'APPRENTISSAGE ET DE VALIDATION
#   Le regulateur n'appelle le PINN qu'avec la charge nominale et la tension
#   d'entree de son observateur. Domaine : R de 4 a 6 ohms (plage nominale
#   +-20 %, celle de la boite des gains), Vin de 150 a 240 V, Vout de 0 a
#   200 V, iL de 0 a 40 A ; donnees : transitions de ELM_A1 et ELM_A2 dans
#   ce domaine ; validation : fenetres de ELM_V1 et ELM_V2 entierement dans
#   ce domaine.
#
# CRITERES
#   Fenetres de 0.5 ms (110 pas Tc), prediction recurrente du PINN depuis
#   l'etat vrai au debut de la fenetre, avec la fraction de conduction et
#   les vrais Vin et G :
#     C1 : sur chaque enregistrement, erreur efficace mediane <= 10 mV et
#          centile 99 <= 20 mV ;
#     C2 : erreur maximale <= 0.5 V sur toutes les fenetres.
#   Candidats : A, 3 couches de 32 neurones (valeur d'Ito et Wasa) ; B,
#   3 couches de 64. Choix : le plus petit qui passe C1 et C2. Si aucun ne
#   passe, on garde le candidat de plus faible erreur mediane comme modele
#   du regulateur, avec son ecart au modele physique.
#   Reference : le modele physique (integre exactement) sur les memes
#   fenetres.
#
# REPRODUCTIBILITE
#   Les tirages sont fixes (donnees : graine 7 ; points de physique :
#   hypercube latin, graine 8 ; poids de depart : Glorot, graine 1).
#   L-BFGS sur plusieurs milliers d'iterations amplifie les ecarts
#   d'arrondi : sur une autre machine (autre bibliotheque de calcul
#   lineaire), les poids peuvent differer dans les derniers chiffres et le
#   chemin d'optimisation peut s'ecarter. Le fichier de reference est
#   pinn_pid_modele.mat livre avec ce script ; une nouvelle execution doit
#   redonner le meme verdict sur les criteres, pas forcement les memes
#   poids.
#
# CE QUE PRODUIT CE SCRIPT
#   pinn_pid_modele.mat : poids du PINN retenu (W1..W4, b1..b4), constantes
#     de normalisation et d'echelle, vecteurs de test (entrees, sorties et
#     derivees d'un pas) que le bloc MATLAB doit retrouver ;
#   entrainement_pinn_resultats.json : reference du modele physique,
#     historique des pertes, criteres de chaque candidat, choix.
#
# FICHIERS NECESSAIRES (dans le meme dossier que ce script)
#   banc_commun.py, scenarios_communs.json, scenario_S*.mat,
#   donnees_elm_ELM_A1.mat, donnees_elm_ELM_A2.mat, donnees_elm_ELM_V1.mat,
#   donnees_elm_ELM_V2.mat, scenario_ELM_A1.mat ... scenario_ELM_V2.mat.
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy)
#
# COMMENT LANCER CE SCRIPT
#   python entrainement_pinn.py
#   Les candidats deja entraines (pinn_candidat_A_nominale.npz et
#   pinn_candidat_B_nominale.npz, livres avec ce script) sont repris
#   s'ils ont les reglages du script : la validation et le choix prennent
#   alors quelques minutes. Sans ces fichiers, l'entrainement prend de
#   deux a quatre heures (B est le plus long). Pour
#   entrainer les deux candidats en parallele sur deux coeurs :
#   "python entrainement_pinn.py A" et "python entrainement_pinn.py B"
#   dans deux consoles, puis "python entrainement_pinn.py".
#   pinn_pid_modele.mat est livre : il n'est pas necessaire de relancer ce
#   script pour la suite.
#
# ORDRE D'EXECUTION (PINN-PID)
#   1. ce script ; 2. estimation_etat_pinn.py ; 3. boite_gains_pinn.py ;
#   4. choix_correction_pinn.py ; 5. banc_pinn_pid.py ;
#   6. Tester_PINN_PID_Rejeu.m ; 7. Construction_PINN_PID.m ;
#   8. Simuler_PINN_PID.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques, plage et definitions du banc commun

import os                                     # chemins de fichiers
import sys                                    # argument de la ligne de commande
from types import SimpleNamespace             # petit objet a attributs
import json                                   # fichier de resultats
import time                                   # duree d'execution
import numpy as np                            # calcul numerique
from scipy.io import loadmat, savemat         # fichiers .mat
from scipy.optimize import minimize           # L-BFGS
from scipy.stats import qmc                   # hypercube latin
from scipy.linalg import expm                 # modele physique de reference

PLAGE = "nominale"                            # domaine d'apprentissage : charge nominale

try:                                          # dossier du script (ou dossier courant)
    DOSSIER_PINN = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_PINN = os.getcwd()

with open(os.path.join(DOSSIER_PINN, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC or "iL_max_dem_A" not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py ne contient pas les definitions attendues (marque de fin, observation du courant).")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])  # constantes du circuit, Circuit, simuler

T_DEBUT = time.time()                         # pour la duree totale
NF = 110                                      # pas Tc par fenetre de 0.5 ms
SV, SI = TC * 10.0 / C_CONV, TC * 100.0 / L_BOB   # echelles de sortie du reseau (V, A)
S = np.array([SV, SI])
N_DONNEES, N_PHYSIQUE = 40000, 40000          # transitions et points de physique
GRAINE_DONNEES, GRAINE_PHYSIQUE, GRAINE_POIDS = 7, 8, 1
ITER_MAX = 15000                              # iterations L-BFGS au plus
CANDIDATS = {"A": 32, "B": 64}                # largeur des 3 couches cachees

R_MIN, R_MAX, IL_MAX = 4.0, 6.0, 40.0         # domaine d'apprentissage
G_CENTRE, G_ECHELLE = 0.2, 0.05               # normalisation de G (centree sur 1/5 S)
VIN_MIN, VIN_MAX, V_MAX = 150.0, 240.0, 200.0


def titre(texte):
    """Bandeau de titre dans la console."""
    print("\n" + "=" * 78 + "\n " + texte + "\n" + "=" * 78)


print(f"Plage {PLAGE} : R de {R_MIN:g} a {R_MAX:g} ohms, Vin de {VIN_MIN:g} a {VIN_MAX:g} V.")


# %% ETAPE 1 : le reseau, ses gradients et la perte

def frac_tranche(d, k):
    """Fraction de la tranche k pendant laquelle le MOSFET conduit (porteuse
    de 1 200 pas, 120 pas par tranche, comme le banc commun)."""
    n_on = np.ceil(d * NPER).astype(int)
    return np.clip(n_on - (k * NREG) % NPER, 0, NREG) / NREG


def normaliser(tau_n, v0, i0, s, vin, G):
    """Entrees du reseau, ramenees vers [-1 ; 1] : tau_n = tau / Tc."""
    return np.stack([2 * tau_n - 1, (v0 - 100.0) / 60.0, (i0 - 12.0) / 10.0, 2 * s - 1,
                     (vin - 195.0) / 45.0, (G - G_CENTRE) / G_ECHELLE], axis=1)


def f_phys(v, i, s, vin, G):
    """Second membre du modele moyen par tranche (avec blocage de la diode)."""
    fv = (i - G * v) / C_CONV
    fi = (s * vin - v - (1 - s) * VF - (s * RON + (1 - s) * RDIODE) * i) / L_BOB
    bloque = (i <= 0) & (fi < 0)                                # la diode empeche le courant de s'inverser
    fi = np.where(bloque, 0.0, fi)
    return fv, fi, bloque


class Reseau:
    """Perceptron a couches tanh ; les parametres sont ranges dans un seul
    vecteur theta (W1, b1, W2, b2, ...), W de taille (entrees, sorties)."""

    def __init__(self, tailles, graine):
        self.tailles = tailles
        rng = np.random.default_rng(graine)
        ps = []
        for a, b in zip(tailles[:-1], tailles[1:]):
            lim = np.sqrt(6.0 / (a + b))                        # initialisation de Glorot
            ps += [rng.uniform(-lim, lim, (a, b)).ravel(), np.zeros(b)]
        self.theta = np.concatenate(ps)

    def deballer(self, th):
        """Liste des couches (W, b) a partir du vecteur theta."""
        out, p = [], 0
        for a, b in zip(self.tailles[:-1], self.tailles[1:]):
            W = th[p:p + a * b].reshape(a, b); p += a * b
            bb = th[p:p + b]; p += b
            out.append((W, bb))
        return out

    def avant(self, th, z, tangente=False):
        """Sortie N(z) (n x 2) et, si tangente, sa derivee par rapport a
        tau_n (mode direct : dz0/dtau_n = 2)."""
        couches = self.deballer(th)
        h = z; hs = [z]; hds = []
        hd = None
        if tangente:
            hd = np.zeros_like(z); hd[:, 0] = 2.0
            hds.append(hd)
        for (W, b) in couches[:-1]:
            h = np.tanh(h @ W + b); hs.append(h)
            if tangente:
                hd = (1 - h * h) * (hd @ W); hds.append(hd)
        W, b = couches[-1]
        N = h @ W + b
        Nd = hd @ W if tangente else None
        return N, Nd, (couches, hs, hds)

    def arriere(self, gN, gNd, cache):
        """Gradient par rapport a theta, a partir des gradients par rapport
        a N et a sa derivee en tau (retropropagation)."""
        couches, hs, hds = cache
        grads = [None] * len(couches)
        W, b = couches[-1]
        h = hs[-1]
        gW = h.T @ gN; gb = gN.sum(0)
        gh = gN @ W.T
        ghd = None
        if gNd is not None:
            hd = hds[-1]
            gW += hd.T @ gNd
            ghd = gNd @ W.T
        grads[-1] = (gW, gb)
        for l in range(len(couches) - 2, -1, -1):
            W, b = couches[l]
            h = hs[l + 1]; hprev = hs[l]
            dsig = 1 - h * h                                    # derivee de tanh
            if ghd is not None:
                hd = hds[l + 1]; hdprev = hds[l]
                ad = hdprev @ W                                 # tangente de la pre-activation
                ga_d = ghd * dsig
                gh = gh + ghd * (-2 * h * ad)                   # a travers la derivee de tanh
                gW_d = hdprev.T @ ga_d
                ghd = ga_d @ W.T
            else:
                gW_d = 0.0
            ga = gh * dsig
            gW = hprev.T @ ga + gW_d
            gb = ga.sum(0)
            gh = ga @ W.T
            grads[l] = (gW, gb)
        return np.concatenate([np.concatenate([g[0].ravel(), g[1]]) for g in grads])


def perte(th, net, Dd, Dp, lam=1.0):
    """Perte L = L_donnees + lam L_physique et son gradient."""
    zd, xd0, xd1 = Dd                                           # donnees : tau_n = 1
    N, _, cache = net.avant(th, zd)
    pred = xd0 + S * N
    rd = (pred - xd1) / S
    Ld = np.mean(np.sum(rd ** 2, 1))
    gN = 2 * rd / len(rd)
    g = net.arriere(gN, None, cache)
    zp, xp0, tau_n, sp, vinp, Gp = Dp                           # physique
    N, Nd, cache = net.avant(th, zp, tangente=True)
    x = xp0 + tau_n[:, None] * S * N
    fv, fi, bl = f_phys(x[:, 0], x[:, 1], sp, vinp, Gp)
    r = N + tau_n[:, None] * Nd - TC * np.stack([fv, fi], 1) / S   # (dx/dt - f) Tc / S
    Lp = np.mean(np.sum(r ** 2, 1))
    gr = 2 * r / len(r) * lam
    dfv_dv = -Gp / C_CONV; dfv_di = 1.0 / C_CONV                # derivees de f par rapport a l'etat
    dfi_dv = np.where(bl, 0.0, -1.0 / L_BOB); dfi_di = np.where(bl, 0.0, -(sp * RON + (1 - sp) * RDIODE) / L_BOB)
    t = tau_n
    gNv = gr[:, 0] - (gr[:, 0] * TC / SV * dfv_dv * t * SV + gr[:, 1] * TC / SI * dfi_dv * t * SV)
    gNi = gr[:, 1] - (gr[:, 0] * TC / SV * dfv_di * t * SI + gr[:, 1] * TC / SI * dfi_di * t * SI)
    gN = np.stack([gNv, gNi], 1)
    gNd = gr * t[:, None]
    g = g + net.arriere(gN, gNd, cache)
    return Ld + lam * Lp, g, Ld, Lp


# %% ETAPE 2 : donnees, points de physique et enregistrements de validation

titre("ETAPE 2 : donnees et points de physique")


def charger(code):
    """Enregistrement Simulink de l'ELM : Vout, iL, d, Vin et G."""
    m = loadmat(os.path.join(DOSSIER_PINN, f"donnees_elm_{code}.mat"))
    sc = loadmat(os.path.join(DOSSIER_PINN, f"scenario_{code}.mat"))
    d = m["d"].ravel().astype(float); v = m["vout"].ravel().astype(float); il = m["il"].ravel().astype(float)
    vin = sc["vin"].ravel(); G = 1 / float(sc["R0"].item()) + sc["gx"].ravel()
    return v, il, d, vin, G


def dans_domaine(vin, G):
    """Instants dont la charge et Vin sont dans le domaine de la plage. Une
    marge relative de 1e-9 garde les bornes exactes : 1/(1/98) vaut
    98.00000000000001 en double precision, et sans cette marge les
    fenetres a 98 ohms tout rond seraient exclues."""
    R = 1 / G
    m = 1e-9
    return ((R >= R_MIN * (1 - m)) & (R <= R_MAX * (1 + m)) &
            (vin >= VIN_MIN * (1 - m)) & (vin <= VIN_MAX * (1 + m)))


def donnees(n, graine):
    """n transitions (etat k, entrees k) -> etat k+1, tirees dans ELM_A1 et ELM_A2."""
    X = []
    for code in ("ELM_A1", "ELM_A2"):
        v, il, d, vin, G = charger(code)
        k = np.arange(len(v) - 1)
        Xc = np.stack([v[:-1], il[:-1], frac_tranche(d[:-1], k), vin[:-1], G[:-1], v[1:], il[1:]], 1)
        X.append(Xc[dans_domaine(vin[:-1], G[:-1])])
    X = np.concatenate(X)
    print(f"  Transitions disponibles dans le domaine : {len(X)} ; tirees : {n} (graine {graine}).")
    rng = np.random.default_rng(graine)
    X = X[rng.choice(len(X), n, replace=False)]
    z = normaliser(np.ones(n), X[:, 0], X[:, 1], X[:, 2], X[:, 3], X[:, 4])
    return z, X[:, 0:2], X[:, 5:7]


def physique(n, graine):
    """n points de physique par hypercube latin dans le domaine. s vaut 0
    (45 %), 1 (45 %) ou une valeur intermediaire (10 %)."""
    u = qmc.LatinHypercube(d=6, seed=graine).random(n)
    tau = u[:, 0]; v0 = V_MAX * u[:, 1]; i0 = IL_MAX * u[:, 2]
    s = np.where(u[:, 3] < 0.45, 0.0, np.where(u[:, 3] > 0.55, 1.0, (u[:, 3] - 0.45) / 0.1))
    vin = VIN_MIN + (VIN_MAX - VIN_MIN) * u[:, 4]
    R = np.exp(np.log(R_MIN) + (np.log(R_MAX) - np.log(R_MIN)) * u[:, 5]); G = 1 / R
    z = normaliser(tau, v0, i0, s, vin, G)
    return z, np.stack([v0, i0], 1), tau, s, vin, G


DD = donnees(N_DONNEES, GRAINE_DONNEES)
DP = physique(N_PHYSIQUE, GRAINE_PHYSIQUE)
print(f"  Points de physique : {N_PHYSIQUE} (hypercube latin, graine {GRAINE_PHYSIQUE}).")

VALIDATION = {}                               # nom -> (v, iL, d, Vin, G)
for code in ("ELM_V1", "ELM_V2"):
    VALIDATION[code] = charger(code)

def fenetres(v, vin, G):
    """Debuts des fenetres de validation : a partir de la 3e, entierement
    dans le domaine de la plage."""
    nfen = (len(v) - 1) // NF
    k0 = np.arange(2, nfen) * NF
    dom = dans_domaine(vin, G)
    return k0[np.array([dom[k:k + NF + 1].all() for k in k0])]


# %% ETAPE 3 : reference, le modele physique sur les fenetres de validation

titre("ETAPE 3 : reference (modele physique integre exactement)")


def pas_physique(v, i, s, G, vin):
    """Un pas Tc du modele moyen par tranche pour plusieurs etats (memes
    formules que estimation_etat_pinn.py)."""
    B = v.shape[0]
    req = s * RON + (1 - s) * RDIODE
    M = np.zeros((B, 3, 3))
    M[:, 0, 0] = -G / C_CONV * TC; M[:, 0, 1] = TC / C_CONV
    M[:, 1, 0] = -TC / L_BOB; M[:, 1, 1] = -req / L_BOB * TC
    M[:, 1, 2] = (s * vin - (1 - s) * VF) / L_BOB * TC
    E = expm(M)
    vn = E[:, 0, 0] * v + E[:, 0, 1] * i + E[:, 0, 2]
    inn = E[:, 1, 0] * v + E[:, 1, 1] * i + E[:, 1, 2]
    neg = (inn < 0) & (s < 1)
    if np.any(neg):
        a = np.clip(np.where(neg, i / np.maximum(i - inn, 1e-12), 1.0), 0.0, 1.0)
        vmid = v + (vn - v) * a
        vn = np.where(neg, vmid * np.exp(-G * (1 - a) * TC / C_CONV), vn)
        inn = np.where(neg, 0.0, inn)
    return vn, inn


def erreurs_fenetres(pas_un, v, il, d, vin, G):
    """Prediction recurrente de NF pas depuis l'etat vrai au debut de chaque
    fenetre de validation ; renvoie l'erreur efficace et l'erreur maximale
    de chaque fenetre (V)."""
    k0 = fenetres(v, vin, G)
    s_tout = frac_tranche(d, np.arange(len(v)))
    x = np.stack([v[k0], il[k0]], 1)
    err = np.zeros((len(k0), NF))
    for j in range(NF):
        kk = k0 + j
        x = pas_un(x, s_tout[kk], vin[kk], G[kk])
        err[:, j] = x[:, 0] - v[kk + 1]
    return np.sqrt(np.mean(err ** 2, 1)), np.abs(err).max(1)


def pas_un_physique(x, s, vin, G):
    vn, inn = pas_physique(x[:, 0], x[:, 1], s, G, vin)
    return np.stack([vn, inn], 1)


REFERENCE = {}
for nom, (v, il, d, vin, G) in VALIDATION.items():
    rms, mx = erreurs_fenetres(pas_un_physique, v, il, d, vin, G)
    REFERENCE[nom] = {"fenetres": int(len(rms)), "med_mV": float(np.median(rms) * 1e3),
                      "p99_mV": float(np.percentile(rms, 99) * 1e3), "max_mV": float(mx.max() * 1e3)}
    print(f"  {nom:7s} {len(rms):4d} fenetres : mediane {REFERENCE[nom]['med_mV']:.2f} mV, centile 99 "
          f"{REFERENCE[nom]['p99_mV']:.2f} mV, pire {REFERENCE[nom]['max_mV']:.1f} mV")


# %% ETAPE 4 : entrainement et validation des candidats

titre("ETAPE 4 : entrainement des candidats (L-BFGS) et criteres C1, C2")


def controle_gradient(net, th):
    """Gradient analytique contre differences finies sur 50 transitions et
    50 points de physique."""
    sub = lambda D, m: tuple(x[:m] for x in D)
    dd, dp = sub(DD, 50), sub(DP, 50)
    _, g0, _, _ = perte(th, net, dd, dp)
    rng = np.random.default_rng(0)
    pire = 0.0
    for idx in rng.choice(len(th), 6, replace=False):
        e = np.zeros_like(th); e[idx] = 1e-6
        fd = (perte(th + e, net, dd, dp)[0] - perte(th - e, net, dd, dp)[0]) / 2e-6
        pire = max(pire, abs(fd - g0[idx]) / max(abs(fd), 1e-8))
    return pire


# Lancement : "python entrainement_pinn.py" entraine les candidats qui
# n'ont pas encore de fichier pinn_candidat_<nom>_<plage>.npz, puis valide et
# choisit. "python entrainement_pinn.py A" entraine seulement le candidat A
# et s'arrete (pour entrainer les deux candidats en parallele sur une
# machine a deux coeurs) ; une execution sans argument reprend ensuite les
# deux fichiers. Un fichier n'est repris que si ses reglages sont ceux de ce
# script (largeur, iterations, graines, plage, normalisation).
ARGS = sys.argv[1:]
SEUL = ARGS[0] if ARGS else None             # candidat a entrainer seul, ou None
if SEUL is not None and SEUL not in CANDIDATS:
    raise ValueError(f"Candidat inconnu : {SEUL} (choix : {', '.join(CANDIDATS)}).")
SIGNATURE = (f"plage={PLAGE};iter={ITER_MAX};graines={GRAINE_DONNEES},{GRAINE_PHYSIQUE},{GRAINE_POIDS};"
             f"n={N_DONNEES},{N_PHYSIQUE};norm=iL-12/10,G-{G_CENTRE}/{G_ECHELLE};iLmax={IL_MAX}")


def entrainer(nom, H):
    """Entraine un candidat par L-BFGS, l'enregistre et renvoie
    (theta, iterations, message, historique des pertes)."""
    net = Reseau([6, H, H, H, 2], GRAINE_POIDS)
    ecart = controle_gradient(net, net.theta)
    print(f"\n  Candidat {nom} ({H} neurones par couche) : controle du gradient, ecart relatif maximal {ecart:.1e}")
    if ecart > 1e-4:
        raise RuntimeError("Le gradient analytique ne correspond pas aux differences finies.")
    t0 = time.time(); hist = []

    def fun(th):
        L, g, Ld, Lp = perte(th, net, DD, DP)
        hist.append((L, Ld, Lp))
        if len(hist) % 500 == 0:
            print(f"    evaluation {len(hist)} : perte {L:.3e} (donnees {Ld:.3e}, physique {Lp:.3e}), "
                  f"{time.time() - t0:.0f} s", flush=True)
        return L, g

    res = minimize(fun, net.theta.copy(), jac=True, method="L-BFGS-B",
                   options={"maxiter": ITER_MAX, "maxfun": int(ITER_MAX * 1.25), "ftol": 0.0, "gtol": 0.0,
                            "maxcor": 20})
    print(f"    fin : {res.message} ; {res.nit} iterations, {time.time() - t0:.0f} s")
    np.savez(os.path.join(DOSSIER_PINN, f"pinn_candidat_{nom}_{PLAGE}.npz"), theta=res.x, nit=res.nit,
             message=str(res.message), hist=np.array(hist), signature=SIGNATURE)
    return res.x, int(res.nit), str(res.message), hist


if SEUL is not None:
    entrainer(SEUL, CANDIDATS[SEUL])
    print(f"  Candidat {SEUL} entraine et enregistre ; relancer sans argument pour la validation et le choix.")
    sys.exit(0)

RESULTATS = {}
POIDS = {}
for nom, H in CANDIDATS.items():
    net = Reseau([6, H, H, H, 2], GRAINE_POIDS)
    fichier = os.path.join(DOSSIER_PINN, f"pinn_candidat_{nom}_{PLAGE}.npz")
    if os.path.exists(fichier) and str(np.load(fichier)["signature"]) == SIGNATURE:
        f_ = np.load(fichier)
        th, nit, message, hist = f_["theta"], int(f_["nit"]), str(f_["message"]), [tuple(h) for h in f_["hist"]]
        print(f"\n  Candidat {nom} ({H} neurones par couche) : repris de {os.path.basename(fichier)} "
              f"({nit} iterations, {message}).")
    else:
        th, nit, message, hist = entrainer(nom, H)
    L, _, Ld, Lp = perte(th, net, DD, DP)
    print(f"    perte {L:.3e} (donnees {Ld:.3e}, physique {Lp:.3e})")
    res = SimpleNamespace(x=th, nit=nit, message=message)   # meme interface que le resultat de minimize
    POIDS[nom] = (H, res.x)
    crit = {}
    for vnom, (v, il, d, vin, G) in VALIDATION.items():
        def pas_un_pinn(x, s, vinj, Gj, th=res.x, net=net):
            z = normaliser(np.ones(len(x)), x[:, 0], x[:, 1], s, vinj, Gj)
            Nn, _, _ = net.avant(th, z)
            return x + S * Nn
        rms, mx = erreurs_fenetres(pas_un_pinn, v, il, d, vin, G)
        c1 = bool(np.median(rms) <= 0.010 and np.percentile(rms, 99) <= 0.020)
        c2 = bool(mx.max() <= 0.5)
        crit[vnom] = {"med_mV": float(np.median(rms) * 1e3), "p99_mV": float(np.percentile(rms, 99) * 1e3),
                      "max_mV": float(mx.max() * 1e3), "C1": c1, "C2": c2}
        print(f"    {vnom:7s} mediane {crit[vnom]['med_mV']:6.2f} mV, centile 99 {crit[vnom]['p99_mV']:7.2f} mV, "
              f"pire {crit[vnom]['max_mV']:8.1f} mV ; C1 {'oui' if c1 else 'non'}, C2 {'oui' if c2 else 'non'}")
    RESULTATS[nom] = {"neurones": H, "iterations": int(res.nit), "message": str(res.message),
                      "perte": float(L), "perte_donnees": float(Ld), "perte_physique": float(Lp),
                      "historique_toutes_500": [list(map(float, h)) for h in hist[::500]],
                      "validation": crit, "passe": all(c["C1"] and c["C2"] for c in crit.values())}


# %% ETAPE 5 : choix et ecriture des fichiers

titre("ETAPE 5 : choix du candidat")
passent = [n for n in CANDIDATS if RESULTATS[n]["passe"]]
if passent:
    CHOIX = passent[0]                        # le plus petit qui passe (ordre de la liste)
    print(f"  Retenu : candidat {CHOIX} (le plus petit qui passe C1 et C2).")
else:
    CHOIX = min(CANDIDATS, key=lambda n: max(c["med_mV"] for c in RESULTATS[n]["validation"].values()))
    print(f"  Aucun candidat ne passe C1 et C2. Garde comme modele : {CHOIX} (plus faible erreur mediane).")
for vnom in VALIDATION:
    c = RESULTATS[CHOIX]["validation"][vnom]; r = REFERENCE[vnom]
    print(f"  {vnom:7s} PINN {c['med_mV']:.2f} / {c['p99_mV']:.2f} / {c['max_mV']:.1f} mV ; modele physique "
          f"{r['med_mV']:.2f} / {r['p99_mV']:.2f} / {r['max_mV']:.1f} mV (mediane / centile 99 / pire)")

sortie = {"plage": PLAGE, "domaine": {"R": [R_MIN, R_MAX], "Vin": [VIN_MIN, VIN_MAX],
          "Vout_max": V_MAX, "iL_max": IL_MAX}, "reference_modele_physique": REFERENCE, "candidats": RESULTATS,
          "choix": CHOIX, "aucun_ne_passe": not passent}
nom_json = "entrainement_pinn_resultats.json"
with open(os.path.join(DOSSIER_PINN, nom_json), "w", encoding="utf-8") as f:
    json.dump(sortie, f, indent=1)
print(f"\n  {nom_json} ecrit.")

if PLAGE == "nominale":
    H, th = POIDS[CHOIX]
    net = Reseau([6, H, H, H, 2], GRAINE_POIDS)
    couches = net.deballer(th)
    # Vecteurs de test : un pas du PINN (tau = Tc) en quelques points, avec
    # les derivees par rapport a Vout, iL et s (mode direct), pour le bloc MATLAB.
    X_TEST = np.array([[100.0, 20.0, 1.0, 200.0, 0.2], [100.0, 20.0, 0.0, 200.0, 0.2],
                       [100.0, 20.0, 0.37, 180.0, 0.2], [50.0, 10.0, 0.5, 230.0, 0.2],
                       [150.0, 30.0, 0.8, 160.0, 0.2], [0.0, 0.0, 1.0, 200.0, 0.2]])
    Y_TEST = np.zeros((len(X_TEST), 2)); D_TEST = np.zeros((len(X_TEST), 3, 2))
    for q, (v0, i0, s0, vin0, G0) in enumerate(X_TEST):
        z = normaliser(np.ones(1), np.array([v0]), np.array([i0]), np.array([s0]), np.array([vin0]), np.array([G0]))[0]
        T = np.zeros((3, 6)); T[0, 1] = 1 / 60.0; T[1, 2] = 1 / 10.0; T[2, 3] = 2.0
        h = z
        for W, b in couches[:-1]:
            h = np.tanh(h @ W + b)
            T = (T @ W) * (1 - h * h)
        W, b = couches[-1]
        N = h @ W + b
        Y_TEST[q] = [v0 + SV * N[0], i0 + SI * N[1]]
        D_TEST[q] = (T @ W) * S[None, :]                        # derivees de [v', iL'] - [v, iL] par rapport a v, iL, s
    m = {f"W{l + 1}": couches[l][0] for l in range(len(couches))}
    m.update({f"b{l + 1}": couches[l][1].reshape(1, -1) for l in range(len(couches))})
    m.update({"NEURONES": float(H), "SV": SV, "SI": SI, "G_CENTRE": G_CENTRE, "G_ECHELLE": G_ECHELLE,
              "TC": TC, "X_TEST": X_TEST, "Y_TEST": Y_TEST, "D_TEST": D_TEST.reshape(len(X_TEST), 6),
              "CANDIDAT": CHOIX, "PLAGE": PLAGE})
    savemat(os.path.join(DOSSIER_PINN, "pinn_pid_modele.mat"), m)
    print("  pinn_pid_modele.mat ecrit (poids, normalisation, vecteurs de test).")
print(f"  Duree totale : {time.time() - T_DEBUT:.0f} s.")
