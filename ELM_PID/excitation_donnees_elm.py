# =============================================================================
# excitation_donnees_elm.py
#
# OBJECTIF
#   Preparer les donnees d'apprentissage et de validation de l'ELM. L'ELM
#   apprend comment la tension de sortie moyenne d'une fenetre de 0.5 ms
#   repond au rapport cyclique moyen de cette fenetre et des deux
#   precedentes (structure de Lu et al. a cinq entrees). Il faut donc
#   enregistrer le convertisseur en boucle ouverte, attaque par un rapport
#   cyclique choisi pour l'identification, sur toute la plage de charge et
#   de tension d'entree des essais communs.
#   Ce script :
#     1. fabrique quatre enregistrements (deux d'apprentissage, deux de
#        validation, tirages differents) : rapport cyclique, Vin et charge ;
#     2. les ecrit au format des essais communs (scenario_ELM_*.mat), pour
#        que charger_scenario.m et le modele Simulink les lisent comme les
#        autres essais, avec en plus le rapport cyclique impose (dexc) ;
#     3. simule chaque enregistrement sur le banc commun et ecrit ce que
#        Simulink doit retrouver (predictions_banc_donnees_elm.mat).
#
# POURQUOI LA BOUCLE OUVERTE
#   Des donnees prises en boucle fermee ne disent presque rien du
#   convertisseur : la commande y depend de la sortie par la loi du
#   regulateur, et une regression retrouve le regulateur au lieu du
#   convertisseur (la regression de donnees prises en boucle fermee rend
#   exactement -1/(P + I Te + D/Te) du PID de generation). Ici le rapport cyclique est impose de l'exterieur.
#
# LA FORME DES ENREGISTREMENTS
#   Fenetre : 0.5 ms = 110 periodes du regulateur Tc = 11 periodes de
#   decoupage. Tous les changements tombent sur un debut de fenetre : sur
#   chaque fenetre, le rapport cyclique, Vin et la charge sont constants.
#   - 0 a 40 ms : demarrage a d = 0.5, Vin = 200 V, R = 5 ohms.
#   - Rapport cyclique : paliers de 1 a 40 fenetres (0.5 a 20 ms), duree
#     tiree uniformement sur son logarithme (autant de paliers courts, qui
#     montrent la reaction rapide, que de paliers longs, qui montrent le
#     regime etabli). La valeur vise une tension V* : 7 fois sur 10 autour
#     de 100 V (loi normale d'ecart type 15 V, bornee a 60-140 V), 3 fois sur
#     10 n'importe ou entre 40 et 160 V (plage des consignes de S9). Le
#     rapport cyclique vaut V*/Vin au debut du palier, borne a [0.05 ; 0.95].
#     Une tension visee ne s'eloigne pas de plus de 30 V de la precedente
#     (40 V pour les tirages larges), pour borner les depassements.
#   - Tension d'entree : segments de 10 a 60 ms ; a chaque segment, une
#     valeur visee uniforme entre 150 et 240 V (plage des essais : S5, S8,
#     S9), atteinte par une rampe de 2 V par fenetre au plus (4 V/ms).
#   - Charge : resistance physique R0 = 98 ohms et charge electronique en
#     parallele (comme Buck_Commun.slx). Segments de 10 a 60 ms ; a chaque
#     segment, une resistance visee : 4 fois sur 10 entre 4 et 6 ohms
#     (nominal +-20 %), sinon tiree uniformement sur son logarithme entre 6
#     et 98 ohms. La conductance y va par une rampe de 0.004 S par fenetre au
#     plus (environ 0.8 A/ms a 100 V).
#   POURQUOI DES RAMPES : en boucle ouverte, un echelon de charge de 4 a
#   98 ohms coupe 24 A dans la bobine ; avec une impedance caracteristique
#   racine(L/C) = 14.6 ohms, la tension monte vers 300 V (essai fait : 315 V
#   avec des echelons). Les rampes gardent la tension entre 40 et 170 V,
#   la plage utile, comme la rampe de charge de S9.
#   - Pas de bruit de mesure, pas de quantification : on identifie le
#     convertisseur, pas le capteur.
#
# LES QUATRE ENREGISTREMENTS
#   code     role           duree   graine
#   ELM_A1   apprentissage  1 s     31
#   ELM_A2   apprentissage  1 s     32
#   ELM_V1   validation     0.5 s   41
#   ELM_V2   validation     0.5 s   42
#   Soit 4 000 fenetres d'apprentissage pour les 18 poids de sortie de
#   l'ELM (12 neurones, 5 liaisons directes, une constante), et 2 000 de
#   validation. 3 s simulees en tout sous Simulink.
#   Les graines sont differentes de celle du bruit des essais (20261002) :
#   aucun enregistrement ne reprend un essai commun.
#
# CE QUE PRODUIT CE SCRIPT
#   scenario_ELM_A1.mat ... scenario_ELM_V2.mat : memes champs que les
#     essais communs (t, dvref, vin, gx, bruit, R0, q, ...) plus dexc,
#     le rapport cyclique impose a chaque instant k*Tc, et graine.
#   predictions_banc_donnees_elm.mat : pour chaque enregistrement, la
#     tension et le courant de la bobine moyens par fenetre de 0.5 ms
#     (v_fen_<code>, iL_fen_<code>), et la tension toutes les millisecondes
#     (v_ms_<code>), calcules par le banc commun.
#   excitation_donnees_elm.png : les profils de ELM_A1 et sa reponse.
#
# BIBLIOTHEQUES NECESSAIRES (pip install numpy scipy matplotlib)
#
# COMMENT LANCER CE SCRIPT :
#   python excitation_donnees_elm.py
#   depuis le dossier de la base commune (il lit banc_commun.py). Une
#   demi-minute.
#
# ORDRE D'EXECUTION (donnees de l'ELM-PID)
#   1. ce script ; 2. Construction_Modele_Donnees_ELM.m ;
#   3. Generer_Donnees_ELM.m.
# =============================================================================


# %% ETAPE 0 : bibliotheques et definitions du banc commun

import os                                     # chemins de fichiers
import numpy as np                            # calcul numerique et tirages
from scipy.io import savemat                  # ecriture des fichiers .mat
import matplotlib.pyplot as plt               # figure

try:                                          # dossier du script
    DOSSIER_ELM = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER_ELM = os.getcwd()

# Le banc commun est execute jusqu'a sa marque de fin des definitions :
# constantes (TC, D_MIN, ...), circuit et boucle fermee (simuler).
with open(os.path.join(DOSSIER_ELM, "banc_commun.py"), encoding="utf-8") as f:
    SOURCE_BANC = f.read()
MARQUE = "# === FIN DES DEFINITIONS DU BANC COMMUN ==="
if MARQUE not in SOURCE_BANC:
    raise RuntimeError("banc_commun.py n'a pas la marque de fin des definitions.")
exec(SOURCE_BANC[:SOURCE_BANC.index(MARQUE)])

NF = 110                                      # instants Tc par fenetre de 0.5 ms
R0_DONNEES = 98.0                             # resistance physique (ohms), la plus grande de la plage
SANS_QUANTIF = 1e-9                           # pas de quantification sans effet (comme les essais)
ENREGISTREMENTS = [("ELM_A1", "apprentissage", 1.0, 31), ("ELM_A2", "apprentissage", 1.0, 32),
                   ("ELM_V1", "validation", 0.5, 41), ("ELM_V2", "validation", 0.5, 42)]


# %% ETAPE 1 : fabrication d'un enregistrement

def paliers(rng, n_fen, duree_min_fen, duree_max_fen, tirer_valeur, log_duree):
    """Suite de paliers sur n_fen fenetres : renvoie la valeur de chaque
    fenetre. La duree de chaque palier (en fenetres) est tiree entre
    duree_min_fen et duree_max_fen, uniformement sur son logarithme si
    log_duree, sinon uniformement ; la valeur par tirer_valeur(rng, debut)."""
    valeurs = np.zeros(n_fen)
    i = 0
    while i < n_fen:
        if log_duree:
            duree = int(round(np.exp(rng.uniform(np.log(duree_min_fen), np.log(duree_max_fen)))))
        else:
            duree = int(rng.integers(duree_min_fen, duree_max_fen + 1))
        valeurs[i:i + duree] = tirer_valeur(rng, i)
        i += duree
    return valeurs


def segments_rampes(rng, n_fen, depart, tirer_cible, pente_max):
    """Valeur par fenetre : segments de 20 a 120 fenetres (10 a 60 ms) ;
    a chaque segment une valeur visee tiree, atteinte depuis la valeur
    courante par une rampe d'au plus pente_max par fenetre."""
    x = np.zeros(n_fen)
    courant, i = depart, 0
    while i < n_fen:
        duree = int(rng.integers(20, 121))   # duree du segment en fenetres
        cible = tirer_cible(rng)
        for j in range(i, min(i + duree, n_fen)):
            courant = courant + float(np.clip(cible - courant, -pente_max, pente_max))
            x[j] = courant
        i += duree
    return x


def fabriquer(code, role, duree, graine):
    """Profils d'un enregistrement, a chaque instant k*Tc (k = 0 ... N-1)."""
    rng = np.random.default_rng(graine)
    n = int(round(duree / TC)) + 1            # instants, 0 compris (comme les essais communs)
    n_fen = (n + NF - 1) // NF                # fenetres couvrant l'enregistrement
    n_dem = int(round(0.040 / TC)) // NF      # fenetres du demarrage (40 ms)
    reste = n_fen - n_dem                     # fenetres apres le demarrage
    # Tension d'entree : valeur visee uniforme 150-240 V, rampe de 2 V par fenetre au plus
    vin_f = segments_rampes(rng, reste, 200.0, lambda r: r.uniform(150.0, 240.0), 2.0)

    # Charge : resistance visee (4 fois sur 10 entre 4 et 6 ohms, sinon log-uniforme 6-98 ohms),
    # conductance atteinte par une rampe de 0.004 S par fenetre au plus
    def conductance_visee(r):
        R = r.uniform(4.0, 6.0) if r.uniform() < 0.4 else np.exp(r.uniform(np.log(6.0), np.log(98.0)))
        return 1.0 / R

    r_f = 1.0 / segments_rampes(rng, reste, 1.0 / 5.0, conductance_visee, 0.004)

    # Rapport cyclique : paliers de 1 a 40 fenetres, tension visee V*, d = V*/Vin au debut du palier
    etat = {"v": 100.0}                       # tension visee du palier precedent

    def valeur_d(r, i):
        if r.uniform() < 0.7:
            v_vise, ecart_max = float(np.clip(r.normal(100.0, 15.0), 60.0, 140.0)), 30.0
        else:
            v_vise, ecart_max = r.uniform(40.0, 160.0), 40.0
        v_vise = float(np.clip(v_vise, etat["v"] - ecart_max, etat["v"] + ecart_max))
        etat["v"] = v_vise
        return float(np.clip(v_vise / vin_f[i], 0.05, 0.95))

    d_f = paliers(rng, reste, 1, 40, valeur_d, True)
    # Demarrage, puis passage des fenetres aux instants Tc
    vin_fen = np.concatenate([np.full(n_dem, 200.0), vin_f])
    r_fen = np.concatenate([np.full(n_dem, 5.0), r_f])
    d_fen = np.concatenate([np.full(n_dem, 0.5), d_f])
    vin = np.repeat(vin_fen, NF)[:n]
    r_eq = np.repeat(r_fen, NF)[:n]
    dexc = np.repeat(d_fen, NF)[:n]
    gx = np.maximum(1.0 / r_eq - 1.0 / R0_DONNEES, 0.0)      # 1/R = 1/R0 + gx
    return {"code": code, "nom": f"Donnees ELM, {role}", "duree": duree, "R0": R0_DONNEES, "q": SANS_QUANTIF,
            "evenements": [], "graine": graine,
            "description": f"boucle ouverte, rapport cyclique impose ; {role} ; graine {graine}",
            "t": np.arange(n) * TC, "dvref": np.zeros(n), "vin": vin, "gx": gx, "bruit": np.zeros(n),
            "dexc": dexc}


class CommandeImposee:
    """Remplace le regulateur dans la boucle du banc : renvoie a l'instant k
    le rapport cyclique impose dexc[k], sans regarder l'erreur."""

    def __init__(self, dexc):
        self.dexc = dexc                      # rapport cyclique a chaque instant k*Tc
        self.k = 0                            # instant courant

    def pas(self, e):
        d = self.dexc[self.k]
        self.k += 1
        return d


def moyennes_fenetres(x):
    """Moyenne de x sur chaque fenetre complete de NF instants."""
    m = len(x) // NF
    return x[:m * NF].reshape(m, NF).mean(axis=1)


# %% ETAPE 2 : ecriture des enregistrements et simulation sur le banc

predictions = {"Te": TC, "NF": float(NF)}
sims = {}
for code, role, duree, graine in ENREGISTREMENTS:
    e = fabriquer(code, role, duree, graine)
    assert np.all(e["gx"] >= 0.0) and np.all((e["dexc"] >= 0.05) & (e["dexc"] <= 0.95)), code
    savemat(os.path.join(DOSSIER_ELM, f"scenario_{code}.mat"),
            {"code": e["code"], "nom": e["nom"], "duree": e["duree"], "Te": TC, "R0": e["R0"], "q": e["q"],
             "evenements": np.zeros((1, 0)), "description": e["description"], "graine": float(graine),
             "t": e["t"].reshape(-1, 1), "dvref": e["dvref"].reshape(-1, 1), "vin": e["vin"].reshape(-1, 1),
             "gx": e["gx"].reshape(-1, 1), "bruit": e["bruit"].reshape(-1, 1), "dexc": e["dexc"].reshape(-1, 1)},
            do_compression=True)
    sim = simuler(CommandeImposee(e["dexc"]), e)                  # boucle ouverte sur le banc commun
    assert np.allclose(sim["d"], e["dexc"]), code                # le banc applique bien d impose
    sims[code] = (e, sim)
    predictions[f"v_fen_{code}"] = moyennes_fenetres(sim["v"]).reshape(-1, 1)
    predictions[f"iL_fen_{code}"] = moyennes_fenetres(sim["iL"]).reshape(-1, 1)
    predictions[f"v_ms_{code}"] = sim["v"][::int(round(1e-3 / TC))].reshape(-1, 1)
    r_eq = 1.0 / (1.0 / e["R0"] + e["gx"])
    print(f"  {code} ({role}, {duree} s, graine {graine}) : Vout de {sim['v'].min():.1f} a {sim['v'].max():.1f} V, "
          f"d de {e['dexc'].min():.3f} a {e['dexc'].max():.3f}, Vin de {e['vin'].min():.0f} a {e['vin'].max():.0f} V, "
          f"R de {r_eq.min():.2f} a {r_eq.max():.1f} ohms, iL < 10 mA sur {np.mean(sim['iL'] < 0.01) * 100:.1f} % "
          f"des instants")
savemat(os.path.join(DOSSIER_ELM, "predictions_banc_donnees_elm.mat"), predictions, do_compression=True)
print("  scenario_ELM_*.mat et predictions_banc_donnees_elm.mat ecrits.")


# %% ETAPE 3 : figure

e, sim = sims["ELM_A1"]
t_ms = e["t"] * 1e3
fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
axes[0].plot(t_ms, e["dexc"], lw=0.6)
axes[0].set_ylabel("d impose")
axes[1].plot(t_ms, e["vin"], lw=0.8)
axes[1].set_ylabel("Vin (V)")
axes[2].semilogy(t_ms, 1.0 / (1.0 / e["R0"] + e["gx"]), lw=0.8)
axes[2].set_ylabel("R (ohms)")
axes[3].plot(t_ms, sim["v"], lw=0.5)
axes[3].set_ylabel("Vout banc (V)")
axes[3].set_xlabel("temps (ms)")
axes[0].set_title("ELM_A1 : rapport cyclique, tension d'entree et charge imposes ; reponse du banc commun",
                  fontsize=9, loc="left")
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(DOSSIER_ELM, "excitation_donnees_elm.png"), dpi=110)
plt.close(fig)
print("  excitation_donnees_elm.png ecrit.")
