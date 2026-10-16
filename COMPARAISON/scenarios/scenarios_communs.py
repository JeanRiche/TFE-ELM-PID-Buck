# =============================================================================
# scenarios_communs.py
#
# OBJECTIF
#   Definir, une fois pour toutes et pour toutes les methodes du memoire, les
#   onze essais de developpement (S1 a S9, dont S7a, S7b, S8a et S8b), et les ecrire dans des fichiers
#   .mat que liront a l'identique Simulink (modele commun Buck_Commun.slx) et
#   le banc Python (banc_commun.py). Les sequences de bruit sont tirees ici,
#   avec une graine fixee : toutes les methodes voient exactement le meme
#   bruit.
#
# LES ESSAIS
#   S1   nominal, R = 5 ohms : reponse dynamique (montee, depassement,
#        etablissement, ondulation).
#   S2   F1 a la consigne : consigne = 100 + 2 exp(cos(10 pi t)) sur
#        [0.05 ; 0.07) s (reproduction de l'article de Lu et al.).
#   S3   F2 : charge de 4.333 ohms sur [0.05 ; 0.07) s (article).
#   S4   charge +-20 % : 4 ohms sur [0.05 ; 0.07) s, 6 ohms sur [0.10 ; 0.12) s.
#   S5   Vin +-20 % : 160 V sur [0.05 ; 0.07) s, 240 V sur [0.10 ; 0.12) s.
#   S6   rampe de consigne : 100 -> 110 V sur [0.05 ; 0.07) s, palier, puis
#        110 -> 100 V sur [0.10 ; 0.12) s.
#   S7a  quantification de la mesure par un convertisseur analogique-
#        numerique de 12 bits (pleine echelle 2.5 V derriere un pont 1/80,
#        soit 200 V a l'echelle de Vout, pour couvrir les 150 V de S9) :
#        pas de 2.5 / 4096 x 80 = 48.8 mV a l'echelle de Vout.
#   S7b  bruit gaussien sur la mesure : 10 mV d'ecart type, filtre par un
#        filtre du premier ordre a 10 kHz (filtre anti-repliement), graine
#        fixee.
#   S8a  charge legere, R = 25 ohms pendant tout l'essai, avec Vin +-20 %.
#   S8b  charge tres legere, R = 98 ohms pendant tout l'essai, avec Vin +-20 %.
#   S9   profil grand signal (dynamique et non-linearites du Buck), 0.42 s :
#          0 a 50 ms      demarrage sous 5 ohms (rapport cyclique en butee) ;
#          50 a 90 ms     rampe de charge de 5 a 98 ohms (le courant de charge
#                         baisse de 20 A a 1 A ; l'amortissement du filtre LC
#                         passe de 1.46 a 0.074) ;
#          120 ms         consigne 100 -> 150 V a charge legere ;
#          150 ms         consigne 150 -> 50 V (rapport cyclique en butee
#                         basse, conduction discontinue : la diode se bloque) ;
#          210 ms         Vin 200 -> 150 V (le gain du convertisseur, d x Vin,
#                         change) ;
#          270 ms         consigne 50 -> 100 V avec Vin = 150 V ;
#          330 a 370 ms   Vin revient a 200 V et rampe de charge de 98 a 5 ohms
#                         (le courant remonte de 1 A a 20 A).
#
# COMMENT LA CHARGE VARIE (meme mecanisme dans Simulink et sur le banc)
#   La charge est une resistance physique R0 (fixe pendant un essai : la plus
#   grande resistance de l'essai), en parallele avec une charge electronique
#   en mode "resistance constante" : une source de courant commandee qui tire
#   i = gx(t) x Vout. gx(t) est une conductance (en siemens), toujours >= 0.
#   Charge equivalente : 1 / R = 1 / R0 + gx. Une charge electronique est le
#   materiel qu'on utilise en laboratoire pour faire varier une charge sans
#   a-coup ; elle permet les rampes de S9, impossibles avec des
#   interrupteurs (un echelon de 5 a 98 ohms envoie l'energie de la bobine
#   dans le condensateur, quel que soit le regulateur).
#
# CE QUE CHAQUE FICHIER CONTIENT
#   scenario_<code>.mat, une variable par grandeur :
#     code, nom, duree (s), Te (s) = periode du regulateur Tc = 1/220000 s,
#     R0 (ohms), q (pas de quantification en V,
#     1e-9 = pas de quantification), evenements (instants en s),
#     t (s), dvref (V, ajoute a la consigne de 100 V), vin (V), gx (S),
#     bruit (V, ajoute a la mesure de Vout).
#   scenarios_communs.json : la liste des essais, leurs descriptions, leurs
#     evenements et les choix (graine, niveau de bruit, pas du CAN).
#   scenarios_communs.png : les profils des onze essais.
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT :
#   python scenarios_communs.py
#   Duree : quelques secondes. Les fichiers produits (scenario_<code>.mat et
#   scenarios_communs.json) sont ceux que lisent toutes les methodes, dans
#   leurs dossiers respectifs.
# =============================================================================


# %% ETAPE 0 : bibliotheques et constantes

import os                                     # chemins de fichiers
import json                                   # liste des essais au format texte
import numpy as np                            # tableaux et tirages aleatoires
from scipy.io import savemat                  # ecriture des fichiers .mat
import matplotlib.pyplot as plt               # figure des profils

TE = 1.0 / 220000.0                           # pas commun des profils = periode du regulateur Tc
                                              # (dix echantillons par periode de decoupage de 22 kHz,
                                              # 120 pas du powergui de 1/(22000*1200) s)
VREF = 100.0                                  # consigne nominale (V)
VIN_NOM = 200.0                               # tension d'entree nominale (V)
R_NOM = 5.0                                   # charge nominale (ohms)
R_F2 = 5.0 + (-2.0 / 15.0) * 5.0              # charge de F2 (article) : 4.333 ohms
GRAINE_BRUIT = 20261002                       # graine du bruit de S7b (fixee une fois)
ECART_TYPE_BRUIT = 0.010                      # 10 mV d'ecart type
FC_BRUIT = 10e3                               # frequence de coupure du filtre anti-repliement (Hz)
PAS_CAN = 2.5 / 4096 * 80.0                   # 12 bits, 2.5 V, pont 1/80 : 48.8 mV (pleine echelle 200 V)
SANS_QUANTIF = 1e-9                           # pas "nul" : la quantification ne change rien

try:                                          # dossier du script
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
print("Dossier de travail :", DOSSIER)


# %% ETAPE 1 : outils pour decrire un profil

def pas_de(instant):
    """Numero du pas correspondant a un instant (en secondes)."""
    return int(round(instant / TE))


def palier(n, valeur_base, paliers):
    """Profil constant par morceaux : valeur_base partout, sauf sur chaque
    fenetre (debut, fin, valeur) de la liste 'paliers'. Les fenetres sont
    definies en numeros de pas entiers, pour que Simulink et le banc
    basculent exactement au meme pas."""
    x = np.full(n, float(valeur_base))
    for debut, fin, valeur in paliers:
        x[pas_de(debut):pas_de(fin)] = valeur
    return x


def rampe(x, debut, fin, v0, v1):
    """Remplace x sur [debut ; fin) par une rampe lineaire de v0 a v1."""
    a, b = pas_de(debut), pas_de(fin)
    x[a:b] = v0 + (v1 - v0) * np.arange(b - a) / (b - a)
    return x


def conductance(r_equivalente, r0):
    """Conductance gx que doit tirer la charge electronique pour que la
    charge totale (R0 en parallele) vaille r_equivalente."""
    return 1.0 / r_equivalente - 1.0 / r0


# %% ETAPE 2 : les onze essais

def essai(code, nom, duree, R0, evenements, description, dvref=None, vin=None, gx=None, bruit=None, q=SANS_QUANTIF):
    """Rassemble un essai. Les profils absents prennent leur valeur nominale."""
    n = pas_de(duree) + 1                                  # nombre de pas, instant 0 compris
    return {"code": code, "nom": nom, "duree": duree, "R0": R0, "q": q,
            "evenements": list(evenements), "description": description,
            "t": np.arange(n) * TE,
            "dvref": np.zeros(n) if dvref is None else dvref(n),
            "vin": np.full(n, VIN_NOM) if vin is None else vin(n),
            "gx": np.zeros(n) if gx is None else gx(n),
            "bruit": np.zeros(n) if bruit is None else bruit(n)}


def F1(n):
    """F1 de Lu et al. sur [0.05 ; 0.07) s, en numeros de pas entiers."""
    x = np.zeros(n)
    k = np.arange(pas_de(0.05), pas_de(0.07))
    x[k] = 2.0 * np.exp(np.cos(10.0 * np.pi * k * TE))
    return x


def bruit_filtre(n):
    """Bruit gaussien de graine fixee, filtre au premier ordre a FC_BRUIT,
    ramene a un ecart type de ECART_TYPE_BRUIT sur toute la duree."""
    blanc = np.random.default_rng(GRAINE_BRUIT).standard_normal(n)
    a = np.exp(-2.0 * np.pi * FC_BRUIT * TE)               # pole du filtre discret
    y = np.zeros(n)
    etat = 0.0
    for k in range(n):                                     # y(k) = a y(k-1) + (1 - a) blanc(k)
        etat = a * etat + (1.0 - a) * blanc[k]
        y[k] = etat
    return ECART_TYPE_BRUIT * y / np.std(y)


def vin_pm20(n):
    """Vin : 160 V sur [0.05 ; 0.07), 240 V sur [0.10 ; 0.12), 200 V ailleurs."""
    return palier(n, VIN_NOM, [(0.05, 0.07, 160.0), (0.10, 0.12, 240.0)])


def dvref_rampe(n):
    """Consigne : rampe de 0 a +10 V sur [0.05 ; 0.07), palier, retour a 0 sur [0.10 ; 0.12)."""
    x = palier(n, 0.0, [(0.07, 0.10, 10.0)])
    rampe(x, 0.05, 0.07, 0.0, 10.0)
    return rampe(x, 0.10, 0.12, 10.0, 0.0)


R0_S9 = 98.0                                               # resistance physique de S9
G5_S9 = conductance(R_NOM, R0_S9)                          # conductance pour 5 ohms equivalents


def gx_S9(n):
    """Charge de S9 : 5 ohms, rampe vers 98 ohms sur [0.05 ; 0.09), 98 ohms,
    rampe vers 5 ohms sur [0.33 ; 0.37), puis 5 ohms. La rampe porte sur la
    conductance : le courant de charge varie lineairement a tension fixe."""
    x = palier(n, 0.0, [(0.0, 0.05, G5_S9), (0.37, 10.0, G5_S9)])
    rampe(x, 0.05, 0.09, G5_S9, 0.0)
    return rampe(x, 0.33, 0.37, 0.0, G5_S9)


ESSAIS = [
    essai("S1", "Nominal", 0.2, R_NOM, [], "R = 5 ohms, aucune perturbation"),
    essai("S2", "F1 a la consigne", 0.2, R_NOM, [0.05, 0.07],
          "consigne = 100 + F1(t), F1 = 2 exp(cos(10 pi t)) sur [0.05 ; 0.07) s", dvref=F1),
    essai("S3", "F2 (article)", 0.2, R_NOM, [0.05, 0.07], "charge de 4.333 ohms sur [0.05 ; 0.07) s",
          gx=lambda n: palier(n, 0.0, [(0.05, 0.07, conductance(R_F2, R_NOM))])),
    essai("S4", "Charge +-20 %", 0.2, 6.0, [0.05, 0.07, 0.10, 0.12],
          "4 ohms sur [0.05 ; 0.07) s, 6 ohms sur [0.10 ; 0.12) s, 5 ohms ailleurs",
          gx=lambda n: palier(n, conductance(R_NOM, 6.0), [(0.05, 0.07, conductance(4.0, 6.0)), (0.10, 0.12, 0.0)])),
    essai("S5", "Vin +-20 %", 0.2, R_NOM, [0.05, 0.07, 0.10, 0.12],
          "Vin = 160 V sur [0.05 ; 0.07) s, 240 V sur [0.10 ; 0.12) s", vin=vin_pm20),
    essai("S6", "Rampe de consigne", 0.2, R_NOM, [0.05, 0.07, 0.10, 0.12],
          "consigne 100 -> 110 V sur [0.05 ; 0.07) s, palier, 110 -> 100 V sur [0.10 ; 0.12) s", dvref=dvref_rampe),
    essai("S7a", "Quantification 12 bits", 0.2, R_NOM, [],
          "mesure quantifiee au pas de 48.8 mV (CAN 12 bits, 2.5 V, pont 1/80)", q=PAS_CAN),
    essai("S7b", "Bruit de mesure", 0.2, R_NOM, [],
          "bruit gaussien de 10 mV filtre a 10 kHz sur la mesure, graine fixee", bruit=bruit_filtre),
    essai("S8a", "Charge legere 25 ohms", 0.2, 25.0, [0.05, 0.07, 0.10, 0.12],
          "R = 25 ohms pendant tout l'essai, Vin = 160 V puis 240 V", vin=vin_pm20),
    essai("S8b", "Charge legere 98 ohms", 0.2, 98.0, [0.05, 0.07, 0.10, 0.12],
          "R = 98 ohms pendant tout l'essai, Vin = 160 V puis 240 V", vin=vin_pm20),
    essai("S9", "Profil grand signal", 0.42, R0_S9, [0.05, 0.09, 0.12, 0.15, 0.21, 0.27, 0.33, 0.37],
          "demarrage 5 ohms ; rampe 5 -> 98 ohms ; consigne 150 V, 50 V ; Vin 150 V ; consigne 100 V ; "
          "Vin 200 V et rampe 98 -> 5 ohms",
          dvref=lambda n: palier(n, 0.0, [(0.12, 0.15, 50.0), (0.15, 0.27, -50.0)]),
          vin=lambda n: palier(n, VIN_NOM, [(0.21, 0.33, 150.0)]), gx=gx_S9),
]
print(f"{len(ESSAIS)} essais definis.")


# %% ETAPE 3 : controles et ecriture des fichiers

for e in ESSAIS:
    r_eq = 1.0 / (1.0 / e["R0"] + e["gx"])                 # charge equivalente a chaque pas
    assert np.all(e["gx"] >= 0.0), e["code"]               # une charge electronique ne fournit pas de courant
    assert np.all(np.isfinite(e["dvref"] + e["vin"] + e["gx"] + e["bruit"])), e["code"]
    texte_q = "aucune" if e["q"] == SANS_QUANTIF else f"pas de {e['q'] * 1e3:.1f} mV"   # pour l'affichage
    print(f"  {e['code']:4s} {e['nom']:24s} {e['duree']*1e3:5.0f} ms, {len(e['t'])} pas ; R equivalente de "
          f"{r_eq.min():.3f} a {r_eq.max():.3f} ohms ; consigne de {VREF + e['dvref'].min():.2f} a "
          f"{VREF + e['dvref'].max():.2f} V ; Vin de {e['vin'].min():.0f} a {e['vin'].max():.0f} V ; "
          f"bruit {np.std(e['bruit'])*1e3:.1f} mV ; quantification : {texte_q}")
    savemat(os.path.join(DOSSIER, f"scenario_{e['code']}.mat"),
            {"code": e["code"], "nom": e["nom"], "duree": e["duree"], "Te": TE, "R0": e["R0"], "q": e["q"],
             "evenements": np.array(e["evenements"], dtype=float).reshape(1, -1),
             "description": e["description"], "t": e["t"].reshape(-1, 1),
             "dvref": e["dvref"].reshape(-1, 1), "vin": e["vin"].reshape(-1, 1),
             "gx": e["gx"].reshape(-1, 1), "bruit": e["bruit"].reshape(-1, 1)},
            do_compression=True)

liste = {"format": 2, "Te": TE, "Ts_powergui": 1.0 / (22000.0 * 1200.0), "Vref": VREF, "Vin_nominale": VIN_NOM, "graine_bruit": GRAINE_BRUIT,
         "ecart_type_bruit_V": ECART_TYPE_BRUIT, "fc_bruit_Hz": FC_BRUIT, "pas_CAN_V": PAS_CAN,
         "bande_retour_V": 1.0,
         "essais": [{k: e[k] for k in ("code", "nom", "duree", "R0", "q", "evenements", "description")} for e in ESSAIS]}
with open(os.path.join(DOSSIER, "scenarios_communs.json"), "w", encoding="utf-8") as f:
    json.dump(liste, f, indent=1)
print("  scenario_*.mat et scenarios_communs.json ecrits.")


# %% ETAPE 4 : figure des profils

fig, axes = plt.subplots(len(ESSAIS), 1, figsize=(11, 2.0 * len(ESSAIS)), sharex=False)
for ax, e in zip(axes, ESSAIS):
    t_ms = e["t"] * 1e3
    ax.plot(t_ms, VREF + e["dvref"], label="consigne (V)")                  # consigne vue par le comparateur
    ax.plot(t_ms, e["vin"] / 2.0, label="Vin / 2 (V)")                      # tension d'entree, divisee par 2
    if np.any(e["bruit"]):                                                   # bruit amplifie 100 fois pour etre visible
        ax.plot(t_ms, 100.0 + 100.0 * e["bruit"], lw=0.3, label="100 + 100 x bruit (V)")
    ax2 = ax.twinx()                                                         # second axe pour la charge
    ax2.semilogy(t_ms, 1.0 / (1.0 / e["R0"] + e["gx"]), "k:", label="R equivalente (ohms, axe droit)")
    ax2.set_ylim(2.0, 200.0)
    ax.set_title(f"{e['code']} : {e['nom']}", fontsize=9, loc="left")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, loc="upper left")
    ax2.legend(fontsize=6, loc="upper right")
axes[-1].set_xlabel("temps (ms)")
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER, "scenarios_communs.png"), dpi=110)
plt.close(fig)
print("  scenarios_communs.png ecrit.")
