# =============================================================================
# schemas_chapitre4.py
#
# OBJET
#   Dessiner les trois schemas du chapitre 4 (conception de l'ELM-PID) :
#     Fig4_1_boucle_ELM_PID        boucle fermee et boucle d'adaptation ;
#     Fig4_2_reseau_ELM            structure du reseau ELM et de son jacobien ;
#     Fig4_3_chronogramme_fenetres pas Tc, fenetres de 0,5 ms, porte, zone morte.
#   Ce sont des schemas : aucune donnee simulee n'y figure. Chaque element
#   dessine correspond a une operation de ELM_PID/banc_elm_pid.py (classes
#   BlocPID, AdaptateurELM, ELMPID), de ELM_PID/elm_pid_adaptatif.m et de
#   ELM_PID/banc_commun.py (fonction simuler) ; les reglages sont ceux de
#   ELM_PID/banc_elm_pid.py (dictionnaire REGLAGES).
#
# COMMENT LANCER
#   python schemas_chapitre4.py
#   Sorties dans le dossier du script : PNG 300 dpi et PDF vectoriel.
#   Prerequis : LaTeX (latex, dvipng) avec les paquets newtx et icomma ; les textes et
#   formules sont composes par LaTeX (text.usetex), dans les polices du memoire.
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Arc

DOSSIER = os.path.dirname(os.path.abspath(__file__))
CM = 1 / 2.54
LARGEUR = 15.0                      # cm, bloc de texte de 15,5 cm
ROUGE = (0.80, 0.10, 0.10)          # couleur de l'ELM-PID au chapitre 5
GRIS = (0.45, 0.45, 0.45)
GRIS_CLAIR = (0.93, 0.93, 0.93)
TAILLE = 10

# Texte et formules composes par LaTeX avec les polices du memoire (Times :
# newtxtext et newtxmath) ; la virgule decimale des formules vient d'icomma.
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{newtxtext,newtxmath}"
                           r"\usepackage{icomma}",
    "font.family": "serif",
    "font.size": TAILLE,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def nouvelle_figure(hauteur):
    fig = plt.figure(figsize=(LARGEUR * CM, hauteur * CM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, LARGEUR)
    ax.set_ylim(0, hauteur)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def boite(ax, x0, y0, x1, y1, texte="", couleur="black", lw=0.9, taille=TAILLE, fond="white", z=3):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="square,pad=0",
                                fc=fond, ec=couleur, lw=lw, zorder=z))
    if texte:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, texte, ha="center", va="center",
                fontsize=taille, zorder=z + 1, linespacing=1.25)


def fleche(ax, pts, couleur="black", lw=0.9, tete=True):
    """Polyligne terminee par une tete de fleche."""
    pts = np.asarray(pts, float)
    if len(pts) > 2:
        ax.plot(pts[:-1, 0], pts[:-1, 1], color=couleur, lw=lw, zorder=2, solid_capstyle="butt")
    if tete:
        ax.annotate("", xy=pts[-1], xytext=pts[-2],
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.18",
                                    color=couleur, lw=lw, shrinkA=0, shrinkB=0), zorder=2)
    else:
        ax.plot(pts[-2:, 0], pts[-2:, 1], color=couleur, lw=lw, zorder=2)


def point(ax, x, y, couleur="black"):
    ax.add_patch(Circle((x, y), 0.055, fc=couleur, ec=couleur, zorder=5))


def sauver(fig, nom):
    for ext, opts in (("png", {"dpi": 300}), ("pdf", {})):
        fig.savefig(os.path.join(DOSSIER, f"{nom}.{ext}"), facecolor="white",
                    metadata={"Creator": None, "Producer": None} if ext == "pdf" else {"Software": None}, **opts)
    plt.close(fig)
    print(f"{nom} : png et pdf ecrits")


# -----------------------------------------------------------------------------
# Figure 4.1 : boucle fermee et boucle d'adaptation
# -----------------------------------------------------------------------------
def figure_4_1():
    H = 9.6
    fig, ax = nouvelle_figure(H)
    ym = 2.45                              # ligne principale
    yr = 0.85                              # ligne de retour
    ysep = 5.6                             # separation des deux cadences

    # deux echelles de temps
    ax.plot([0.1, LARGEUR - 0.1], [ysep, ysep], color=GRIS, lw=0.7, ls=(0, (5, 3)), zorder=0)
    ax.text(0.15, H - 0.12, "En fin de chaque fenêtre $n$ de 0,5 ms (110 $T_c$)", ha="left", va="top")
    ax.text(0.15, 0.1, "À chaque pas $T_c$ = 1/220 000 s", ha="left", va="bottom")

    # consigne et sommateur
    ax.text(0.15, ym + 0.12, r"$v_\mathrm{ref}$", ha="left", va="bottom")
    fleche(ax, [(0.15, ym), (0.98, ym)])
    ax.add_patch(Circle((1.2, ym), 0.22, fc="white", ec="black", lw=0.9, zorder=3))
    ax.text(0.86, ym + 0.32, "+", ha="center", va="center", fontsize=9)
    ax.text(0.95, ym - 0.45, "\u2212", ha="center", va="center", fontsize=9)

    # e -> PID
    fleche(ax, [(1.42, ym), (2.0, ym)])
    ax.text(1.68, ym + 0.08, "$e$", ha="left", va="bottom")
    point(ax, 1.62, ym)
    boite(ax, 2.0, ym - 0.7, 4.3, ym + 0.7, "Bloc PID\nparallèle\n($P$, $I$, $D$)")

    # saturation
    fleche(ax, [(4.3, ym), (4.7, ym)])
    boite(ax, 4.7, ym - 0.5, 5.8, ym + 0.5)
    ax.plot([4.85, 5.07, 5.43, 5.65], [ym - 0.3, ym - 0.3, ym + 0.3, ym + 0.3], color="black", lw=0.9, zorder=4)
    ax.text(5.25, ym - 0.58, "[0,01 ; 0,99]", ha="center", va="top")

    # d -> MLI -> Buck
    fleche(ax, [(5.8, ym), (6.5, ym)])
    point(ax, 6.2, ym)
    ax.text(6.1, ym + 0.08, "$d$", ha="right", va="bottom")
    boite(ax, 6.5, ym - 0.5, 7.75, ym + 0.5, "MLI\n22 kHz")
    fleche(ax, [(7.75, ym), (8.1, ym)])
    boite(ax, 8.1, ym - 0.7, 10.5, ym + 0.7, "Convertisseur\nBuck")
    fleche(ax, [(9.3, yr + 0.35), (9.3, ym - 0.7)], couleur=GRIS)
    ax.text(9.42, yr + 0.42, r"$V_\mathrm{in}$, charge", ha="left", va="center")

    # v_o -> mesure -> y
    fleche(ax, [(10.5, ym), (11.1, ym)])
    ax.text(10.8, ym + 0.08, "$v_o$", ha="center", va="bottom")
    boite(ax, 11.1, ym - 0.55, 13.0, ym + 0.55, "Mesure\n(bruit, CAN)")
    fleche(ax, [(13.0, ym), (14.7, ym)])
    point(ax, 13.5, ym)
    ax.text(14.3, ym + 0.08, "$y$", ha="center", va="bottom")
    fleche(ax, [(13.5, ym), (13.5, yr), (1.2, yr), (1.2, ym - 0.22)])

    # moyennes par fenetre (cumul a chaque Tc, sortie en fin de fenetre)
    ya0, ya1 = 3.55, 5.0
    boite(ax, 5.0, ya0, 13.0, ya1,
          "Moyennes sur la fenêtre $n$ : $\\bar{e}(n)$, $\\bar{y}(n)$, $\\bar{d}(n)$\n"
          "et sensibilité $\\sigma(n)$ de $d$ aux gains\n"
          "fenêtre suspecte si $d$ en butée à un pas de la fenêtre", couleur=ROUGE)
    xg = 2.75                                # verticale des gains
    ye = 4.25
    ax.plot([1.62, 1.62, xg - 0.15], [ym, ye, ye], color=ROUGE, lw=0.9, zorder=2)
    ax.add_patch(Arc((xg, ye), 0.3, 0.3, theta1=0, theta2=180, color=ROUGE, lw=0.9, zorder=2))
    fleche(ax, [(xg + 0.15, ye), (5.0, ye)], couleur=ROUGE)
    fleche(ax, [(6.2, ym), (6.2, ya0)], couleur=ROUGE)
    fleche(ax, [(13.5, ym), (13.5, ye), (13.0, ye)], couleur=ROUGE)

    # chaine de fin de fenetre (de droite a gauche)
    yb0, yb1 = 6.15, 8.6
    boite(ax, 11.6, yb0, 14.75, yb1,
          "Réseau ELM\nprédiction $\\hat{y}(n)$\njacobien $J(n)$\nOS-ELM si\n$|\\bar{y}-\\hat{y}|$ > 0,1 V", couleur=ROUGE)
    boite(ax, 7.55, yb0, 11.1, yb1,
          "Porte : aucune butée\nen $n$, $n-1$, $n-2$\nZone morte :\n$|\\bar{e}|$ > 0,1 V et $J$ > 0",
          couleur=ROUGE)
    boite(ax, 3.45, yb0, 7.05, yb1,
          "Gradient normalisé\n$\\varphi = J\\,(\\sigma \\odot K_\\mathrm{ZN})$\n"
          "$\\Delta \\theta = \\eta\\,\\bar{e}\\,\\varphi/(\\varepsilon+\\varphi\\cdot\\varphi)$\n"
          "$+\\;\\alpha\\,(\\theta-\\theta_\\mathrm{préc})$", couleur=ROUGE)
    boite(ax, 0.3, yb0, 3.05, yb1, "Projection sur\nl'ensemble\nadmissible\n(Rosen)", couleur=ROUGE)
    # liaisons de la chaine
    yh, yl = yb0 + 1.75, yb0 + 0.55
    fleche(ax, [(11.6, yh), (11.1, yh)], couleur=ROUGE)
    fleche(ax, [(11.1, yl), (11.6, yl)], couleur=ROUGE, lw=0.7)
    fleche(ax, [(7.55, (yb0 + yb1) / 2), (7.05, (yb0 + yb1) / 2)], couleur=ROUGE)
    fleche(ax, [(3.45, (yb0 + yb1) / 2), (3.05, (yb0 + yb1) / 2)], couleur=ROUGE)
    # entrees venant des moyennes
    for xv, lab in ((12.4, "$\\bar{y}$, $\\bar{d}$"), (9.3, "$\\bar{e}$, butée"), (5.25, "$\\bar{e}$, $\\sigma$")):
        fleche(ax, [(xv, ya1), (xv, yb0)], couleur=ROUGE)
        ax.text(xv + 0.1, 5.3, lab, ha="left", va="center")
    # gains vers le PID
    fleche(ax, [(xg, yb0), (xg, ym + 0.7)], couleur=ROUGE)
    ax.text(xg - 0.12, 5.3, "$K = K_\\mathrm{ZN} \\odot \\theta$", ha="right", va="center")
    sauver(fig, "Fig4_1_boucle_ELM_PID")


# -----------------------------------------------------------------------------
# Figure 4.2 : reseau ELM (RVFL : liaisons directes)
# -----------------------------------------------------------------------------
def figure_4_2():
    H = 10.4
    fig, ax = nouvelle_figure(H)
    xe, xc, xs = 2.2, 5.3, 8.3
    noms = ["$\\bar{y}(n-1)$", "$\\bar{y}(n-2)$", "$\\bar{d}(n)$", "$\\bar{d}(n-1)$", "$\\bar{d}(n-2)$"]
    ye = np.linspace(7.0, 2.6, 5)
    # deux groupes de 6 : les liaisons directes passent entre eux
    yc = np.concatenate([np.linspace(9.0, 6.4, 6), np.linspace(3.2, 0.6, 6)])
    ys = 4.8
    r_e, r_c, r_s = 0.22, 0.2, 0.3

    for a in ye:                                   # poids d'entree, fixes
        for b in yc:
            ax.plot([xe + r_e, xc - r_c], [a, b], color=(0.75, 0.75, 0.75), lw=0.35, zorder=1)
    for b in yc:                                   # poids de sortie, appris
        ax.plot([xc + r_c, xs - r_s], [b, ys], color=ROUGE, lw=0.6, zorder=1)
    for a in ye:                                   # liaisons directes, apprises
        ax.plot([xe + r_e, xs - r_s], [a, ys], color=ROUGE, lw=0.8, ls=(0, (4, 2.5)), zorder=2)

    for a, nom in zip(ye, noms):
        ax.add_patch(Circle((xe, a), r_e, fc="white", ec="black", lw=0.9, zorder=3))
        ax.text(xe - 0.32, a, nom, ha="right", va="center")
    for b in yc:
        ax.add_patch(Circle((xc, b), r_c, fc="white", ec="black", lw=0.9, zorder=3))
        sx = np.linspace(-0.12, 0.12, 30)
        ax.plot(xc + sx, b + 0.22 * (1 / (1 + np.exp(-sx * 40)) - 0.5), color="black", lw=0.6, zorder=4)
    ax.add_patch(Circle((xs, ys), r_s, fc="white", ec="black", lw=0.9, zorder=3))
    ax.text(xs, ys, "$\\Sigma$", ha="center", va="center", zorder=4)
    fleche(ax, [(xs + r_s, ys), (9.15, ys)])
    ax.text(9.25, ys, "$\\hat{y}(n)$", ha="left", va="center")

    ax.add_patch(Circle((xs, 7.3), 0.2, fc="white", ec="black", lw=0.9, zorder=3))
    ax.text(xs, 7.3, "1", ha="center", va="center", zorder=4)
    fleche(ax, [(xs, 7.1), (xs, ys + r_s)], couleur=ROUGE, lw=0.6)
    ax.text(xs + 0.12, 6.4, "$\\beta_0$", ha="left", va="center")

    # titres de colonnes
    ax.text(xe - 0.6, H - 0.15, "Entrées\n(centrées réduites)", ha="center", va="top")
    ax.text(xc, H - 0.15, "12 neurones\nlogistiques $g_j$", ha="center", va="top")
    ax.text(xs, H - 0.15, "Sortie\nlinéaire", ha="center", va="top")

    # etiquettes des poids
    ax.text(xe - 0.4, 1.9, "$w_j$, $b_j$ : tirés\nau hasard dans\n[−1 ; 1], fixés",
            ha="center", va="top")
    ax.text(7.0, 7.9, "$\\beta_j$", ha="center", va="center")
    ax.text(xc, 4.8, "$\\beta_{L,i}$", ha="center", va="center", color="black",
            bbox=dict(fc="white", ec="none", pad=0.5), zorder=4)
    ax.text(10.45, 3.55, "$\\hat{y}(n)$ : prédiction de $\\bar{y}(n)$", ha="left", va="center")
    # legende des traits
    lx = 10.45
    ax.plot([lx, lx + 0.7], [2.55, 2.55], color=(0.75, 0.75, 0.75), lw=0.8)
    ax.text(lx + 0.85, 2.55, "poids fixes", ha="left", va="center")
    ax.plot([lx, lx + 0.7], [2.0, 2.0], color=ROUGE, lw=0.8)
    ax.text(lx + 0.85, 2.0, "poids appris", ha="left", va="center")
    ax.plot([lx, lx + 0.7], [1.45, 1.45], color=ROUGE, lw=0.8, ls=(0, (4, 2.5)))
    ax.text(lx + 0.85, 1.45, "liaison directe (apprise)", ha="left", va="center")

    # encadres : apprentissage et jacobien
    ax.text(lx, 9.6,
            "Poids appris $\\beta$ :\n"
            "moindres carrés régularisés\n"
            "hors ligne ($C_\\mathrm{r}$ = 0,1), puis\n"
            "mise à jour OS-ELM en ligne",
            ha="left", va="top")
    ax.text(lx, 7.3,
            "Jacobien :\n"
            "$J(n) = \\partial\\hat{y}(n)/\\partial\\bar{d}(n)$,\n"
            "dérivée exacte du réseau\n"
            "(termes $g_j(1-g_j)$ et\n"
            "liaison directe de $\\bar{d}(n)$)",
            ha="left", va="top")
    sauver(fig, "Fig4_2_reseau_ELM")


# -----------------------------------------------------------------------------
# Figure 4.3 : chronogramme des fenetres
# -----------------------------------------------------------------------------
def figure_4_3():
    H = 11.0
    fig = plt.figure(figsize=(LARGEUR * CM, H * CM))
    gauche, droite = 0.19, 0.975
    nf = 6                                     # fenetres m-1 ... m+4
    noms = ["$m-1$", "$m$", "$m+1$", "$m+2$", "$m+3$", "$m+4$"]
    pas = 11                                   # un trait pour 10 Tc : 11 par fenetre
    t = (np.arange(nf * pas) + 0.5) / pas * 0.5   # ms
    lignes = [("Pas du PID", 0.8), ("Erreur", 1.5), ("Rapport\ncyclique", 1.3),
              ("Fenêtre\nsuspecte", 0.5), ("Porte", 0.95), ("Décision en fin\nde fenêtre", 0.75)]
    tot = sum(h for _, h in lignes)
    bas, haut = 0.12, 0.99
    axes = []
    y = haut
    for nom, h in lignes:
        hh = (haut - bas) * h / tot
        y -= hh
        a = fig.add_axes([gauche, y + 0.01, droite - gauche, hh - 0.02])
        a.set_xlim(0, nf * 0.5 + 0.06)
        for c in ("top", "right", "left"):
            a.spines[c].set_visible(False)
        a.spines["bottom"].set_color(GRIS)
        a.spines["bottom"].set_bounds(0, nf * 0.5)
        a.set_yticks([])
        a.set_xticks([])
        a.set_ylim(0, 1)
        for k in range(nf + 1):
            a.axvline(0.5 * k, color=GRIS, lw=0.5, ls=(0, (3, 3)), zorder=0)
        fig.text(gauche - 0.015, y + hh / 2, nom, ha="right", va="center", linespacing=1.15)
        axes.append(a)
    a_tc, a_e, a_d, a_s, a_p, a_g = axes

    # 1. pas du PID et noms des fenetres
    a_tc.vlines(t, 0, 0.35, color="black", lw=0.6)
    for i in range(nf):
        a_tc.text(0.5 * i + 0.25, 0.5, "fenêtre " + noms[i], ha="center", va="bottom", fontsize=9)

    # 2. erreur e(k) et moyenne ebar(n) ; bande de la zone morte
    ebar = [0.05, 1.20, 0.45, 0.28, 0.18, 0.06]
    e = []
    for i in range(nf):
        tt = np.linspace(0, 1, pas)
        if i == 1:
            v = 2.2 * np.sin(np.pi * np.clip(tt * 1.15, 0, 1))
            v = v - v.mean() + ebar[i]
        else:
            v = ebar[i] + 0.08 * np.sin(2 * np.pi * tt * 1.5 + i)
        e.extend(v)
    e = np.array(e)
    a_e.set_ylim(-0.5, 2.6)
    a_e.axhspan(-0.1, 0.1, xmax=3.0 / 3.06, color=GRIS_CLAIR, zorder=0)
    a_e.plot(t, e, "o", ms=2.2, color="black", zorder=3)
    for i in range(nf):
        a_e.plot([0.5 * i + 0.02, 0.5 * (i + 1) - 0.02], [ebar[i]] * 2, color="black", lw=1.6, zorder=4)
    a_e.text(3.0, 2.5, "points : $e(k)$ ; traits : $\\bar{e}(n)$", ha="right", va="top", fontsize=9,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=5)
    a_e.annotate("zone morte\n$|\\bar{e}|$ < 0,1 V", xy=(0.25, 0.08), xytext=(0.25, 0.75),
                 ha="center", va="bottom", fontsize=9, linespacing=1.1,
                 arrowprops=dict(arrowstyle="-", color=GRIS, lw=0.6))

    # 3. rapport cyclique, butee dans la fenetre m
    d = 0.5 + 0.015 * np.sin(np.arange(nf * pas) * 1.3)
    i0, i1 = pas + 1, pas + 6
    d[i0:i1] = 0.99
    d[i1:i1 + 3] = [0.8, 0.65, 0.56]
    a_d.set_ylim(0.3, 1.12)
    a_d.axhline(0.99, xmax=3.0 / 3.06, color=GRIS, lw=0.6, ls=(0, (4, 2)))
    a_d.text(0.03, 1.0, "0,99", ha="left", va="bottom", fontsize=9)
    a_d.step(t, d, where="mid", color="black", lw=0.9)
    a_d.annotate("$d$ en butée", xy=(t[i0 + 2], 0.99), xytext=(1.25, 1.07), ha="left", va="center",
                 fontsize=9, bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=5,
                 arrowprops=dict(arrowstyle="-", color=GRIS, lw=0.6))

    # 4. fenetre suspecte
    a_s.add_patch(Rectangle((0.5, 0.15), 0.5, 0.6, fc=GRIS, ec="black", lw=0.6))
    a_s.text(1.03, 0.45, "fenêtre $m$ suspecte", ha="left", va="center", fontsize=9,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=5)

    # 5. porte, lue en fin de fenetre : fermee pour m, m+1, m+2
    ouverte = [True, False, False, False, True, True]
    for i in range(nf):
        if ouverte[i]:
            a_p.text(0.5 * i + 0.25, 0.25, "ouverte", ha="center", va="center", fontsize=9)
        else:
            a_p.add_patch(Rectangle((0.5 * i + 0.02, 0.05), 0.46, 0.4, fc="white", ec="black", lw=0.6,
                                    hatch="////"))
            a_p.text(0.5 * i + 0.25, 0.25, "fermée", ha="center", va="center", fontsize=9,
                     bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none"), zorder=5)
    a_p.annotate("", xy=(2.5, 0.62), xytext=(1.0, 0.62),
                 arrowprops=dict(arrowstyle="<->,head_length=0.35,head_width=0.15", color="black", lw=0.7))
    a_p.text(1.75, 0.68, "réouverture 1,5 ms après la fin de $m$", ha="center", va="bottom", fontsize=9,
             bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=5)

    # 6. decision en fin de fenetre
    textes = ["zone\nmorte", "porte\nfermée", "porte\nfermée", "porte\nfermée", "gains\nmodifiés", "zone\nmorte"]
    for i in range(nf):
        xf = 0.5 * (i + 1)
        if i == 4:
            a_g.plot(xf, 0.8, "o", ms=6, color=ROUGE, zorder=5, clip_on=False)
        else:
            a_g.plot(xf, 0.8, "o", ms=6, mfc="white", mec="black", mew=0.8, zorder=5, clip_on=False)
        a_g.text(0.5 * i + 0.25, 0.42, textes[i], ha="center", va="center", fontsize=9, linespacing=1.0)

    # axe du temps
    a_g.spines["bottom"].set_color("black")
    a_g.set_xticks(0.5 * np.arange(nf + 1))
    a_g.set_xticklabels([f"{0.5 * k:.1f}".replace(".", ",") for k in range(nf + 1)])
    a_g.tick_params(axis="x", length=3, width=0.6)
    a_g.set_xlabel("Temps (ms)", labelpad=2)
    sauver(fig, "Fig4_3_chronogramme_fenetres")


if __name__ == "__main__":
    figure_4_1()
    figure_4_2()
    figure_4_3()
