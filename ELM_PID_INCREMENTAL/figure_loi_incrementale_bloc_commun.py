# =============================================================================
# figure_loi_incrementale_bloc_commun.py
#
# OBJECTIF
#   Figure du chapitre 4 (section « Pourquoi le bloc PID commun », tableau
#   « Loi incrementale de Lu et bloc PID commun ») : tension de sortie du
#   convertisseur sur les trois essais du tableau, pour trois regulateurs :
#     - le PID de Ziegler-Nichols (point de depart des deux ELM-PID) ;
#     - l'ELM-PID a loi PID incrementale de Lu (dossier ELM_PID_INCREMENTAL) ;
#     - l'ELM-PID a bloc PID commun (dossier ELM_PID).
#   Essais : (a) S1, demarrage ; (b) S2, saut de consigne ; (c) S4, echelon de
#   charge. Ce sont ceux du tableau. Boucle d'adaptation et reglages identiques
#   pour les deux ELM-PID : seule la loi de commande differe.
#
# PRINCIPE
#   Les simulations reprennent exactement celles de
#   COMPARAISON/loi_incrementale_contre_bloc_commun.py, qui a produit le tableau : le
#   banc de chaque dossier (banc_elm_pid.py, jusqu'a l'etape 2) est execute
#   dans un espace de noms propre, puis chaque regulateur est simule sur
#   l'essai avec la fonction simuler du banc commun. Aucun parametre n'est
#   modifie.
#
# CONTROLE AVANT LA FIGURE
#   Les grandeurs du tableau sont recalculees et comparees aux valeurs
#   imprimees dans le memoire. Si un ecart depasse l'arrondi affiche
#   (demi-unite du dernier chiffre), le script s'arrete avant de tracer.
#
# FENETRES ET GRANDEURS TRACEES (toutes dans les fenetres de mesure du tableau)
#   (a) S1 : tension de sortie de 0 a 6 ms, dans la fenetre de mesure 0-30 ms
#       (IAE de demarrage) ; le depassement est mesure avant 50 ms, mais les
#       pics ont lieu a 1,8 et 1,9 ms, et la tension est revenue a la consigne
#       avant 6 ms.
#   (b) S2 : tension de sortie de 45 a 80 ms, autour du saut de consigne de
#       50 ms, dans la fenetre de mesure 30-200 ms.
#   (c) S4 : IAE cumulee, integrale de |e| de 30 ms a t, de 30 a 200 ms (toute
#       la fenetre de mesure). Les valeurs finales sont les IAE apres 30 ms du
#       tableau ; elles sont verifiees avant le trace et ecrites a droite.
#
# STYLES
#   Ceux de COMPARAISON/figures_chapitre5/Figures_Chapitre5.m : Ziegler-
#   Nichols gris [0,65 0,65 0,65], 2,6 pt ; ELM-PID rouge [0,80 0,10 0,10],
#   0,9 pt ; consigne pointillee grise [0,5 0,5 0,5], 0,8 pt. La loi
#   incrementale prend une couleur distincte (violet), 1,6 pt. Trace de la
#   plus epaisse a la plus fine, pour qu'une courbe confondue avec une autre
#   reste visible.
#
# POLICES
#   LaTeX (text.usetex), avec le preambule du memoire : T1, newtxtext,
#   newtxmath, icomma. Graduations 9 pt, axes 10 pt, virgule decimale.
#
# COMMENT LANCER CE SCRIPT
#   python figure_loi_incrementale_bloc_commun.py
#   (une dizaine de secondes ; le dossier ELM_PID a cote de ELM_PID_INCREMENTAL ;
#   necessite latex, matplotlib, numpy, scipy.) Sortie :
#   figures/fig4_05_loi_incrementale_bloc_commun.pdf (dans ce dossier)
# =============================================================================

import os
import sys
import json
import contextlib
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
RACINE = os.path.dirname(DOSSIER)
SORTIE = os.path.join(DOSSIER, "figures", "fig4_05_loi_incrementale_bloc_commun.pdf")

# Noms des regulateurs : cle -> (dossier du banc ou None, etiquette de legende)
ZN, LU, BLOC = "Ziegler-Nichols", "Loi incrementale", "Bloc commun"
ETIQUETTES = {ZN: "Ziegler-Nichols",
              LU: "ELM-PID, loi incr\\'ementale",
              BLOC: "ELM-PID, bloc commun"}

# Styles : couleur et epaisseur (pt) ; Ziegler-Nichols et ELM-PID comme dans
# Figures_Chapitre5.m, loi incrementale en violet (distincte du rouge, y
# compris en vision des couleurs deuteranope).
STYLES = {ZN: ((0.65, 0.65, 0.65), 2.6),
          LU: ((0.35, 0.15, 0.70), 1.6),
          BLOC: ((0.80, 0.10, 0.10), 0.9)}
GRIS_CONSIGNE = (0.5, 0.5, 0.5)
ORDRE_TRACE = (ZN, LU, BLOC)             # la plus epaisse d'abord

# Essais : code, lettre du panneau, grandeur tracee ("v" : tension de sortie ;
# "iae" : IAE cumulee depuis 30 ms), fenetre (ms), limites et graduations de
# l'axe vertical, graduations du temps (ms)
ESSAIS = (("S1", "(a)", "v", (0.0, 6.0), (92.0, 116.0), (95, 100, 105, 110, 115), (0, 2, 4, 6)),
          ("S2", "(b)", "v", (45.0, 80.0), (96.0, 104.0), (96, 98, 100, 102, 104), (50, 60, 70, 80)),
          ("S4", "(c)", "iae", (30.0, 200.0), (0.0, 27.0), (0, 10, 20), (30, 80, 130, 180)))

# Valeurs imprimees dans le tableau du chapitre 4 : (essai, grandeur, nombre
# de decimales, ZN, loi incrementale, bloc commun)
TABLEAU = (("S4", "IAE apres 30 ms (mV.s)", 2, 23.61, 19.57, 18.87),
           ("S2", "IAE apres 30 ms (mV.s)", 2, 2.84, 5.61, 2.41),
           ("S2", "ecart maximal (V)", 2, 2.01, 5.34, 2.01),
           ("S1", "IAE de 0 a 30 ms (mV.s)", 2, 93.77, 81.12, 93.72),
           ("S1", "depassement (%)", 1, 14.0, 2.4, 14.0))


# -----------------------------------------------------------------------------
# 1. Simulations (memes fonctions que COMPARAISON/loi_incrementale_contre_bloc_commun.py)
# -----------------------------------------------------------------------------

def charger(dossier):
    """Definitions du banc d'un dossier (jusqu'a l'etape 2 de banc_elm_pid.py),
    dans un espace de noms propre."""
    chemin = os.path.join(RACINE, dossier, "banc_elm_pid.py")
    ns = {"__file__": chemin, "__name__": "definitions"}
    with open(chemin, encoding="utf-8") as f:
        src = f.read()
    with contextlib.redirect_stdout(io.StringIO()):
        exec(src[:src.index("# %% ETAPE 2 :")], ns)
    return ns


def lancer(ns, fabrique, code):
    """Simule un regulateur sur un essai ; renvoie la simulation et les
    grandeurs du banc, avec l'IAE de 0 a 30 ms."""
    sc = {s["code"]: s for s in ns["SCENARIOS"]}[code]
    sim = ns["simuler"](fabrique(), sc)
    o = ns["grandeurs"](sim, sc)
    k30 = int(round(0.03 / ns["TC"]))
    e = sim["consigne"] - sim["v"]
    return {"sim": sim, "o": o, "IAE_dem": float(np.sum(np.abs(e[:k30])) * ns["TC"]), "TC": ns["TC"]}


def calculer():
    """Les trois regulateurs sur S1, S2 et S4."""
    ns_lu = charger("ELM_PID_INCREMENTAL")
    ns_bloc = charger("ELM_PID")
    res = {}
    for code, *_ in ESSAIS:
        res[code] = {ZN: lancer(ns_bloc, lambda: ns_bloc["PIDClassique"](), code),
                     LU: lancer(ns_lu, lambda: ns_lu["ELMPID"](), code),
                     BLOC: lancer(ns_bloc, lambda: ns_bloc["ELMPID"](), code)}
        # Les deux ELM-PID doivent redonner l'IAE de leur fichier de resultats
        for cle, dossier in ((LU, "ELM_PID_INCREMENTAL"), (BLOC, "ELM_PID")):
            with open(os.path.join(RACINE, dossier, "banc_elm_pid_resultats.json"), encoding="utf-8") as f:
                attendu = json.load(f)["references"]["ELM-PID"][code]["grandeurs"]["IAE"]
            if abs(res[code][cle]["o"]["IAE"] - attendu) > 1e-12:
                raise RuntimeError(f"{cle}, {code} : IAE {res[code][cle]['o']['IAE']} au lieu de {attendu}.")
    return res


def grandeur(r, nom):
    """Valeur d'une grandeur du tableau, dans l'unite imprimee."""
    if nom.startswith("IAE apres 30 ms"):
        return r["o"]["IAE"] * 1e3
    if nom.startswith("IAE de 0 a 30 ms"):
        return r["IAE_dem"] * 1e3
    if nom.startswith("ecart maximal"):
        return r["o"]["e_max_V"]
    if nom.startswith("depassement"):
        return r["o"]["depassement_pct"]
    raise KeyError(nom)


def verifier(res):
    """Compare les valeurs recalculees a celles du tableau ; renvoie True si
    chaque ecart tient dans l'arrondi affiche (demi-unite du dernier chiffre)."""
    ok = True
    print(f"{'Essai':5s} {'Grandeur':26s} {'Regulateur':17s} {'recalcule':>11s} {'imprime':>9s} {'ecart':>9s}")
    for code, nom, nd, *imprimes in TABLEAU:
        for cle, imp in zip((ZN, LU, BLOC), imprimes):
            val = grandeur(res[code][cle], nom)
            ecart = val - imp
            dans = abs(ecart) <= 0.5 * 10.0 ** (-nd) + 1e-9
            ok &= dans
            print(f"{code:5s} {nom:26s} {cle:17s} {val:11.4f} {imp:9.{nd}f} {ecart:+9.4f}  "
                  f"{'ok' if dans else 'ECART > ARRONDI'}")
    return ok


# -----------------------------------------------------------------------------
# 2. Figure
# -----------------------------------------------------------------------------

PREAMBULE = (r"\usepackage[T1]{fontenc}"
             r"\usepackage{newtxtext,newtxmath}"
             r"\usepackage{icomma}")


def virgule(valeur, _=None):
    """Graduation avec virgule decimale, en mode mathematique (icomma)."""
    texte = f"{valeur:g}".replace(".", "{,}")
    return "$" + texte + "$"


def iae_cumulee(r):
    """IAE cumulee de 30 ms a t (mV.s), aux instants du banc : somme de |e| Tc
    sur les instants k >= k30, comme grandeurs() du banc commun."""
    sim = r["sim"]
    k30 = int(round(0.03 / r["TC"]))
    e = sim["consigne"] - sim["v"]
    return sim["t"][k30:] * 1e3, np.cumsum(np.abs(e[k30:])) * r["TC"] * 1e3


def verifier_cumulee(res):
    """Les valeurs finales de l'IAE cumulee de S4 doivent etre celles du
    tableau (IAE apres 30 ms), a l'arrondi pres."""
    ok = True
    imprimes = {ZN: TABLEAU[0][3], LU: TABLEAU[0][4], BLOC: TABLEAU[0][5]}
    assert TABLEAU[0][0] == "S4"
    for cle in ORDRE_TRACE:
        finale = iae_cumulee(res["S4"][cle])[1][-1]
        banc = res["S4"][cle]["o"]["IAE"] * 1e3
        dans = abs(finale - imprimes[cle]) <= 0.5e-2 + 1e-9 and abs(finale - banc) <= 1e-9
        ok &= dans
        print(f"S4 IAE cumulee finale, {cle:17s} {finale:9.4f}  banc {banc:9.4f}  imprime {imprimes[cle]:.2f}  "
              f"{'ok' if dans else 'ECART'}")
    return ok


def figure(res):
    plt.rcParams.update({
        "text.usetex": True,
        "text.latex.preamble": PREAMBULE,
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 3,
    })
    largeur, hauteur = 16.0 / 2.54, 6.0 / 2.54
    fig, axes = plt.subplots(1, 3, figsize=(largeur, hauteur))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.20, top=0.80, wspace=0.42)
    traits = {}                                  # courbes du panneau (a), reprises par la legende

    for ax, (code, lettre, trace, (t0, t1), ylim, yticks, xticks) in zip(axes, ESSAIS):
        r_ref = res[code][ZN]
        t = r_ref["sim"]["t"] * 1e3
        k = (t >= t0 - 1e-9) & (t <= t1 + 1e-9)
        if trace == "v":
            trait = ax.plot(t[k], r_ref["sim"]["consigne"][k], ":", color=GRIS_CONSIGNE, lw=0.8, zorder=1)[0]
            traits.setdefault("consigne", trait)
        for cle in ORDRE_TRACE:
            couleur, epaisseur = STYLES[cle]
            if trace == "v":
                x, y = t[k], res[code][cle]["sim"]["v"][k]
            else:
                x, y = iae_cumulee(res[code][cle])
            trait = ax.plot(x, y, "-", color=couleur, lw=epaisseur, zorder=2)[0]
            traits.setdefault(cle, trait)
        if trace == "iae":
            # Valeurs finales, en petit, a droite des courbes : ZN et loi
            # incrementale au-dessus de leur courbe, bloc commun au-dessous
            # (les deux dernieres valeurs different de 0,7 mV.s seulement).
            for cle, haut, couleur_texte in ((ZN, True, (0.40, 0.40, 0.40)), (LU, True, STYLES[LU][0]),
                                             (BLOC, False, STYLES[BLOC][0])):
                finale = iae_cumulee(res[code][cle])[1][-1]
                ax.annotate("$" + f"{finale:.2f}".replace(".", "{,}") + "$", (t1 - 1.0, finale),
                            xytext=(0, 3.5 if haut else -3.5), textcoords="offset points",
                            ha="right", va="bottom" if haut else "top", fontsize=8, color=couleur_texte)
        ax.set_xlim(t0, t1)
        ax.set_ylim(*ylim)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_title(lettre, fontsize=9, pad=3)
        ax.set_xlabel("temps (ms)", labelpad=1.5)
        ax.set_ylabel(r"$v_{\mathrm{o}}$ (V)" if trace == "v" else r"IAE cumul\'ee (mV$\cdot$s)", labelpad=1.5)
        ax.xaxis.set_major_formatter(FuncFormatter(virgule))
        ax.yaxis.set_major_formatter(FuncFormatter(virgule))
        ax.grid(True, lw=0.3, alpha=0.5)
        ax.set_axisbelow(True)

    # Legende commune : poignees reprises des courbes tracees (meme couleur,
    # meme epaisseur, meme style de trait).
    cles = ORDRE_TRACE + ("consigne",)
    noms = [ETIQUETTES[c] for c in ORDRE_TRACE] + ["consigne"]
    fig.legend([traits[c] for c in cles], noms, loc="upper center", ncol=4, frameon=False, handlelength=2.4,
               columnspacing=1.3, handletextpad=0.5, bbox_to_anchor=(0.5, 1.0), borderaxespad=0.1)
    return fig


def principal():
    res = calculer()
    print("\nControle : valeurs recalculees contre valeurs imprimees dans le tableau du chapitre 4\n")
    ok = verifier(res)
    print()
    ok &= verifier_cumulee(res)
    if not ok:
        print("\nUn ecart depasse l'arrondi affiche : arret, pas de figure.")
        sys.exit(1)
    print("\nTous les ecarts tiennent dans l'arrondi affiche.")
    fig = figure(res)
    os.makedirs(os.path.dirname(SORTIE), exist_ok=True)
    fig.savefig(SORTIE, format="pdf", metadata={"CreationDate": None})
    plt.close(fig)
    print("Ecrit :", SORTIE)


if __name__ == "__main__":
    principal()
