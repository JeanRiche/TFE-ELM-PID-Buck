# =============================================================================
# schema_chapitre3.py
#
# OBJET
#   Dessiner la figure 3.1 du chapitre 3 : le convertisseur Buck de la base
#   commune (Fig3_1_convertisseur_Buck), avec ses grandeurs et la commande
#   du MOSFET par MLI a partir du rapport cyclique d.
#   C'est un schema : aucune donnee simulee n'y figure. Les valeurs ecrites
#   sont celles de ELM_PID/banc_commun.py (etape 0) et du dimensionnement
#   du chapitre 3.
#
# COMMENT LANCER
#   python schema_chapitre3.py      (Python 3, matplotlib, schemdraw)
#   Sorties dans le dossier du script : PNG 300 dpi et PDF vectoriel
#   (polices incorporees), largeur 15 cm, meme style que les schemas du
#   chapitre 4. Prerequis : LaTeX (latex, dvipng) avec les paquets newtx et
#   icomma ; les textes et formules sont composes par LaTeX (text.usetex),
#   dans les polices du memoire.
# =============================================================================

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import schemdraw
import schemdraw.elements as elm

DOSSIER = os.path.dirname(os.path.abspath(__file__))
CM = 1 / 2.54
LARGEUR = 15.0
TAILLE = 10

# Texte et formules composes par LaTeX avec les polices du memoire (Times :
# newtxtext et newtxmath) ; la virgule decimale des formules vient d'icomma.
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{newtxtext,newtxmath}"
                           r"\usepackage{icomma}\DeclareUnicodeCharacter{03A9}{\ensuremath{\Omega}}"
                           r"\DeclareUnicodeCharacter{00B5}{{\fontencoding{TS1}\selectfont\textmu}}",
    "font.family": "serif",
    "font.size": TAILLE,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
schemdraw.config(font="serif", fontsize=TAILLE, lw=0.9, unit=2.4)

with schemdraw.Drawing(show=False) as d:
    H = 4.0                                   # hauteur de la branche haute
    xs, xn, xl, xc, xr, xt = 0.0, 4.2, 7.6, 8.8, 12.0, 13.4
    # Source d'entree, interrupteur (MOSFET), noeud de commutation
    d.add(elm.SourceV().endpoints((xs, 0), (xs, H))
          .label("$V_\\mathrm{in}$", loc="bottom", ofst=0.3))
    d.add(elm.Line().endpoints((xs, H), (1.2, H)))
    d.add(elm.Switch().endpoints((1.2, H), (3.2, H))
          .label("MOSFET\n$R_\\mathrm{on}$ = 0,1 Ω", loc="bottom", ofst=0.25))
    d.add(elm.Line().endpoints((3.2, H), (xn, H)))
    d.add(elm.Dot().at((xn, H)))
    # Diode de roue libre, cathode vers le noeud de commutation
    d.add(elm.Diode().endpoints((xn, 0), (xn, H))
          .label("$V_\\mathrm{f}$ = 0,8 V", loc="bottom", ofst=0.3))
    # Bobine et courant
    d.add(elm.Inductor2(loops=3).endpoints((xn, H), (xl, H))
          .label("$L$ = 10 mH", loc="top", ofst=0.15))
    d.add(elm.Line().endpoints((xl, H), (xt, H)))
    d.add(elm.Arrow(headwidth=0.18, headlength=0.3).endpoints((5.3, H - 0.55), (6.5, H - 0.55))
          .label("$i_L$", loc="bottom", ofst=0.1))
    # Condensateur et charge
    d.add(elm.Dot().at((xc, H)))
    d.add(elm.Capacitor().endpoints((xc, H), (xc, 0))
          .label("$C$ = 47 µF", loc="bottom", ofst=0.3))
    d.add(elm.Dot().at((xr, H)))
    d.add(elm.Resistor().endpoints((xr, H), (xr, 0))
          .label("$R$", loc="bottom", ofst=0.3))
    # Masse et bornes de sortie
    d.add(elm.Line().endpoints((xs, 0), (xt, 0)))
    for x in (xn, xc, xr):
        d.add(elm.Dot().at((x, 0)))
    d.add(elm.Dot(open=True).at((xt, H)))
    d.add(elm.Dot(open=True).at((xt, 0)))
    d.add(elm.Label().at((xt, H - 0.45)).label("+", halign="center"))
    d.add(elm.Label().at((xt, 0.45)).label("−", halign="center"))
    d.add(elm.Label().at((xt + 0.35, H / 2)).label("$v_\\mathrm{o}$", halign="left"))
    # Commande : rapport cyclique d -> MLI -> grille
    for a, b in (((0.9, H + 1.3), (3.5, H + 1.3)), ((3.5, H + 1.3), (3.5, H + 2.2)),
                 ((3.5, H + 2.2), (0.9, H + 2.2)), ((0.9, H + 2.2), (0.9, H + 1.3))):
        d.add(elm.Line().endpoints(a, b))
    d.add(elm.Label().at((2.2, H + 1.75)).label("MLI 22 kHz", halign="center"))
    d.add(elm.Arrow(headwidth=0.18, headlength=0.3).endpoints((2.2, H + 1.3), (2.2, H + 0.35)))
    d.add(elm.Arrow(headwidth=0.18, headlength=0.3).endpoints((-0.9, H + 1.75), (0.9, H + 1.75))
          .label("$d$", loc="top", ofst=0.1))
    fig = d.draw(show=False)

mfig = fig.fig
w, h = mfig.get_size_inches()
mfig.set_size_inches(LARGEUR * CM, LARGEUR * CM * h / w)
for ext, kw in (("png", {"dpi": 300, "metadata": {"Software": None}}), ("pdf", {"metadata": {"Creator": None, "Producer": None}})):
    mfig.savefig(os.path.join(DOSSIER, f"Fig3_1_convertisseur_Buck.{ext}"),
                 bbox_inches="tight", pad_inches=0.05, facecolor="white", **kw)
print("Fig3_1_convertisseur_Buck : PNG et PDF ecrits.")
