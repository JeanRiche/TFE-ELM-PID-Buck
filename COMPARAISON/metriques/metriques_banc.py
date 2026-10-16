# =============================================================================
# metriques_banc.py
#
# OBJECTIF
#   Calculer, sur le banc commun, les metriques descriptives des cinq
#   methodes (Ziegler-Nichols, PSO-PID, Fuzzy-PID, ELM-PID, PINN-PID) sur
#   les cinq scenarios du memoire (S1, S2, S3, S8a, S10). Elles ne changent
#   aucun classement : celui-ci repose sur l'IAE de la fenetre de
#   classement.
#
# NOTATIONS
#   Te = Tc = 1/220000 s (periode du regulateur et des enregistrements).
#   k = 0, 1, ..., N-1 : numero de l'instant k Te (N = 44 001 pour 0.2 s).
#   v(k) : tension de sortie VRAIE (pas la mesure bruitee ou quantifiee).
#   consigne(k) = 100 + dvref(k) (V) ; e(k) = consigne(k) - v(k) (V).
#   d(k) : rapport cyclique qui entre dans le PWM, dans [0.01 ; 0.99].
#   iL(k) : courant de la bobine (A).
#   Evenements : instants de scenario_<code>.mat, k_j = round(t_j / Te).
#   Une fenetre [a ; b) contient les instants k tels que a <= k < b.
#
# METRIQUES
#   1. Demarrage (memes definitions que banc_commun.py, fonction grandeurs) :
#      k_dem = min(round(0.05/Te), k du premier evenement) ;
#      1.1 temps de montee 10-90 % ; 1.2 depassement
#      max(0, max v - 100) / 100 en % ; 1.3 temps d'etablissement dans
#      +-1 V (dernier instant ou |e| > 1 V) ; 1.4 courant crete de la bobine.
#   2. Fenetre de classement [k_a ; k_b) : S1 [0 ; 6600), S2, S3, S8a
#      [6600 ; N), S10 [22000 ; N). tau(k) = (k - k_a) Te.
#      2.1 IAE = somme |e| Te (mV.s) ; 2.2 ISE = somme e^2 Te ; 2.3 ITAE =
#      somme tau |e| Te (compte depuis le debut de la fenetre, a la
#      difference de l'ITAE du banc commun) ; 2.4 ecart maximal ; 2.5 ecart
#      efficace.
#   3. Erreur en regime permanent, sur [max(k_(j-1), k_j - 1100) ; k_j) avant
#      chaque evenement j (5 ms) et sur les 10 dernieres ms (2 200 instants) :
#      3.1 erreur moyenne signee ; 3.2 erreur moyenne absolue = |moyenne de
#      e| (un biais statique, pas la moyenne de |e|) ; 3.3 ondulation de la
#      tension (crete a crete).
#   4. Par evenement, sur [k_j ; k_(j+1)) : 4.1 IAE ; 4.2 ecart maximal ;
#      4.3 temps de retour dans +-1 V, "pas revenu" si |e| depasse encore
#      1 V sur les 10 dernieres ms de la fenetre.
#   5. Commande et conduction, sur la fenetre de classement : 5.1 moyenne de
#      |d(k) - d(k-1)| ; 5.2 part des instants en butee (d <= 0.01 ou
#      d >= 0.99) ; 5.3 part des instants en conduction discontinue
#      (iL < 0.01 A).
#   6. Cout de calcul du regulateur seul : temps de chaque appel reg.pas
#      (time.perf_counter_ns avant et apres l'appel), par methode et par
#      scenario, moyenne et maximum ; pour l'ELM-PID et le PINN-PID, appels
#      ordinaires, appels de fin de fenetre et appels ou l'adaptation des
#      gains a tourne. Ordre de grandeur relatif entre methodes (Python
#      interprete, pas un materiel embarque) ; le maximum inclut les aleas
#      du systeme.
#   Classement affiche : IAE de la fenetre de classement par ordre croissant,
#   "=" entre voisines dont l'IAE differe de moins de 5 % ; S10 : quatre
#   methodes classees, Ziegler-Nichols en reference.
#
# CONTROLES (le script s'arrete si l'un echoue)
#   - IAE de classement identique a celle des fichiers de resultats
#     (COMPARAISON/banc_*_resultats.json : IAE_dem pour S1, grandeurs.IAE
#     pour S2, S3, S8a ; COMPARAISON/S10/resultats_S10.json : IAE de 100 ms
#     a la fin). Ecart attendu : nul.
#   - Demarrage, evenements et (S2, S3, S8a) commande identiques a la
#     fonction grandeurs du banc (1e-12 relatif).
#   - Traduction Python ligne a ligne de Metriques_Simulink.m
#     (traduction_metriques_simulink.py) appliquee aux memes signaux :
#     memes metriques a 1e-9 pres (relatif, ou absolu sous 1).
#
# DOSSIERS NECESSAIRES
#   ../../PSO_PID, ../../FUZZY_PID, ../../ELM_PID, ../../PINN_PID (complets),
#   ../banc_*_resultats.json, ../S10/resultats_S10.json.
#
# FICHIERS ECRITS (dans ce dossier)
#   metriques_banc.csv (une ligne par metrique, separateur ';'),
#   metriques_banc.md (tableaux prets a copier).
#
# COMMENT LANCER CE SCRIPT
#   python metriques_banc.py
#   (un seul processus, pour que les temps d'execution ne soient pas
#   perturbes par d'autres calculs ; quelques minutes)
# =============================================================================

import os
import sys
import json
import math
import time
import numpy as np

try:
    DOSSIER = os.path.dirname(os.path.abspath(__file__))
except NameError:
    DOSSIER = os.getcwd()
RACINE = os.path.dirname(os.path.dirname(DOSSIER))
COMP = os.path.dirname(DOSSIER)
sys.path.insert(0, DOSSIER)
import traduction_metriques_simulink as TR           # noqa: E402
T_DEBUT = time.time()


def charger(dossier, fichier, marque):
    """Execute les definitions d'un banc (jusqu'a la marque) dans un espace
    de noms propre a la methode ."""
    chemin = os.path.join(RACINE, dossier, fichier)
    ns = {"__file__": chemin, "__name__": "definitions"}
    with open(chemin, encoding="utf-8") as f:
        src = f.read()
    exec(src[:src.index(marque)], ns)
    return ns


NS = {"Ziegler-Nichols": charger("ELM_PID", "banc_elm_pid.py", "# %% ETAPE 2 :"),
      "PSO-PID": charger("PSO_PID", "recherche_pso_pid.py", "\n# %% ETAPE 3 :"),
      "Fuzzy-PID": charger("FUZZY_PID", "banc_fuzzy_pid.py", "# %% ETAPE 2 :"),
      "ELM-PID": None, "PINN-PID": charger("PINN_PID", "banc_pinn_pid.py", "# %% ETAPE 2 :")}
NS["ELM-PID"] = NS["Ziegler-Nichols"]
with open(os.path.join(RACINE, "PSO_PID", "predictions_banc_pso_pid.json"), encoding="utf-8") as f:
    G_PSO = json.load(f)["gains"]
FABRIQUES = {"Ziegler-Nichols": lambda ns: ns["PIDClassique"](),
             "PSO-PID": lambda ns: ns["PIDParallele"](G_PSO["P"], G_PSO["I"], G_PSO["D"]),
             "Fuzzy-PID": lambda ns: ns["FuzzyPID"](),
             "ELM-PID": lambda ns: ns["ELMPID"](),
             "PINN-PID": lambda ns: ns["PINNPID"]()}
METHODES = list(FABRIQUES)
SCENARIOS = ["S1", "S2", "S3", "S8a", "S10"]
NSB = NS["ELM-PID"]
TC = NSB["TC"]
TE_MATLAB = 1 / (22000 * 10)                  # Te de Metriques_Simulink.m (pour la traduction)
for m_ in METHODES:
    if NS[m_]["TC"] != TC or NS[m_]["L_BOB"] != NSB["L_BOB"] or NS[m_]["C_CONV"] != NSB["C_CONV"]:
        raise RuntimeError(f"{m_} : banc commun different (Tc, L ou C).")


def titre(t):
    print("\n" + "=" * 100 + "\n " + t + "\n" + "=" * 100, flush=True)


# %% ETAPE 1 : valeurs attendues (fichiers de resultats)

def lire_refs(nom):
    with open(os.path.join(COMP, nom), encoding="utf-8") as f:
        return json.load(f)["references"]


REFS = {"PSO-PID": lire_refs("banc_pso_pid_resultats.json"), "Fuzzy-PID": lire_refs("banc_fuzzy_pid_resultats.json"),
        "ELM-PID": lire_refs("banc_elm_pid_resultats.json"), "PINN-PID": lire_refs("banc_pinn_pid_resultats.json")}
with open(os.path.join(COMP, "S10", "resultats_S10.json"), encoding="utf-8") as f:
    R_S10 = json.load(f)
ATTENDU = {}                                  # IAE de classement attendue (V.s)
for m_ in METHODES:
    ref = REFS["PSO-PID"]["Ziegler-Nichols"] if m_ == "Ziegler-Nichols" else REFS[m_][m_]
    ATTENDU[m_] = {"S1": ref["S1"]["IAE_dem"], "S2": ref["S2"]["grandeurs"]["IAE"],
                   "S3": ref["S3"]["grandeurs"]["IAE"], "S8a": ref["S8a"]["grandeurs"]["IAE"],
                   "S10_mVs": R_S10["IAE_100_fin_mVs"][m_]["nominal"]}


# %% ETAPE 2 : metriques (numpy, indices a partir de 0 ; ecrit independamment de Metriques_Simulink.m)

K30 = int(round(0.03 / TC))
K100 = int(round(0.10 / TC))
N5 = int(round(0.005 / TC))
N10 = int(round(0.010 / TC))


def fenetre_classement(code, N):
    return {"S1": (0, K30), "S10": (K100, N)}.get(code, (K30, N))


def metriques(sim, sc, code):
    """Liste de (groupe, fenetre, metrique, valeur, unite) ."""
    v, d, iL = sim["v"], sim["d"], sim["iL"]
    e = sim["consigne"] - v
    N = len(v)
    evts = [int(round(x / TC)) for x in sc["evenements"]]
    out = []
    # 1. demarrage
    k_dem = min([int(round(0.05 / TC))] + evts)
    k10, k90 = np.where(v >= 10.0)[0], np.where(v >= 90.0)[0]
    out.append(("demarrage", "demarrage", "t_montee", float((k90[0] - k10[0]) * TC * 1e3) if len(k10) and len(k90)
                else float("nan"), "ms"))
    out.append(("demarrage", "demarrage", "depassement", float(max(0.0, (v[:k_dem].max() - 100.0) / 100.0 * 100.0)), "%"))
    hors = np.where(np.abs(e[:k_dem]) > 1.0)[0]
    out.append(("demarrage", "demarrage", "t_etablissement", float((hors[-1] + 1) * TC * 1e3) if len(hors) else 0.0, "ms"))
    out.append(("demarrage", "demarrage", "iL_crete", float(iL[:k_dem].max()), "A"))
    # 2. fenetre de classement
    ka, kb = fenetre_classement(code, N)
    ew = e[ka:kb]
    tau = np.arange(kb - ka) * TC
    out.append(("classement", "classement", "IAE", float(np.sum(np.abs(ew)) * TC) * 1e3, "mV.s"))
    out.append(("classement", "classement", "ISE", float(np.sum(ew ** 2) * TC) * 1e3, "V2.ms"))
    out.append(("classement", "classement", "ITAE", float(np.sum(tau * np.abs(ew)) * TC) * 1e3, "mV.s2"))
    out.append(("classement", "classement", "e_max", float(np.max(np.abs(ew))), "V"))
    out.append(("classement", "classement", "e_eff", float(np.sqrt(np.mean(ew ** 2))) * 1e3, "mV"))
    # 4. evenements
    ev = []
    for j, kj in enumerate(evts):
        kf = evts[j + 1] if j + 1 < len(evts) else N
        ef = e[kj:kf]
        h = np.where(np.abs(ef) > 1.0)[0]
        revenu = bool(np.max(np.abs(e[max(kj, kf - N10):kf])) <= 1.0)
        t_ret = 0.0 if len(h) == 0 else (float((h[-1] + 1) * TC * 1e3) if revenu else float("nan"))
        ev.append((kj, float(np.sum(np.abs(ef)) * TC) * 1e3, float(np.max(np.abs(ef))), t_ret, revenu))
    n_non = sum(1 for kj, _, _, _, rev in ev if ka <= kj < kb and not rev)
    out.append(("classement", "classement", "n_non_revenus", float(n_non), "-"))
    # 3. regime permanent
    for j, kj in enumerate(evts):
        a = max(evts[j - 1] if j else 0, kj - N5)
        lab = "avant_%gms" % (kj * TC * 1e3)
        em = float(np.mean(e[a:kj]))
        out += [("regime_permanent", lab, "e_moy", em * 1e3, "mV"),
                ("regime_permanent", lab, "e_moy_abs", abs(em) * 1e3, "mV"),
                ("regime_permanent", lab, "ondulation_Vout", float(np.ptp(v[a:kj])) * 1e3, "mV")]
    em = float(np.mean(e[N - N10:]))
    out += [("regime_permanent", "fin", "e_moy", em * 1e3, "mV"),
            ("regime_permanent", "fin", "e_moy_abs", abs(em) * 1e3, "mV"),
            ("regime_permanent", "fin", "ondulation_Vout", float(np.ptp(v[N - N10:])) * 1e3, "mV")]
    for kj, iae, emax, t_ret, rev in ev:
        lab = "ev_%gms" % (kj * TC * 1e3)
        out += [("evenement", lab, "IAE", iae, "mV.s"), ("evenement", lab, "e_max", emax, "V"),
                ("evenement", lab, "t_retour", t_ret, "ms"), ("evenement", lab, "revenu", float(rev), "-")]
    # 5. commande et conduction
    dw = d[ka:kb]
    out.append(("commande", "classement", "dd_moyen", float(np.mean(np.abs(np.diff(dw)))), "-"))
    out.append(("commande", "classement", "butee", float(np.mean((dw <= 0.01 + 1e-12) | (dw >= 0.99 - 1e-12)) * 100.0), "%"))
    out.append(("commande", "classement", "dcm", float(np.mean(iL[ka:kb] < 0.01) * 100.0), "%"))
    return out


# %% ETAPE 3 : chronometrage du regulateur

class Chrono:
    """Enveloppe transparente : appelle reg.pas avec les memes arguments et
    mesure la duree de l'appel. Type d'appel : 0 ordinaire, 1 fin de fenetre
    sans adaptation des gains, 2 fin de fenetre avec adaptation (ligne du
    journal marquee "adapte")."""

    def __init__(self, reg):
        self.reg = reg
        self.utilise_mesure = getattr(reg, "utilise_mesure", False)
        self.jr = getattr(reg, "journal", None)
        self.t_ns, self.type = [], []

    def pas(self, *args):
        n0 = len(self.jr) if self.jr is not None else 0
        t0 = time.perf_counter_ns()
        u = self.reg.pas(*args)
        t1 = time.perf_counter_ns()
        self.t_ns.append(t1 - t0)
        if self.jr is not None and len(self.jr) > n0:
            self.type.append(2 if self.jr[-1].get("adapte", False) else 1)
        else:
            self.type.append(0)
        return u


def stats_temps(t_ns, typ, fenetres):
    t = np.asarray(t_ns, dtype=float) / 1e3   # us
    typ = np.asarray(typ)
    o = {"appel_moy": (t.mean(), "us"), "appel_max": (t.max(), "us"), "n_appels": (float(len(t)), "-")}
    if fenetres:
        for nom, masque in (("ordinaire", typ == 0), ("fin_fenetre", typ >= 1), ("adaptation", typ == 2)):
            n = int(masque.sum())
            o[f"{nom}_moy"] = (float(t[masque].mean()) if n else float("nan"), "us")
            o[f"{nom}_max"] = (float(t[masque].max()) if n else float("nan"), "us")
            o[f"n_{nom}"] = (float(n), "-")
    return o


# %% ETAPE 4 : simulations, controles

titre("ETAPE 4 : simulations (cinq methodes x cinq scenarios, un seul processus) et controles")
BASE = {}                                     # cle -> valeur
LIGNES = []                                   # (scenario, methode, groupe, fenetre, metrique, valeur, unite)
TEMPS = {m_: {"t": [], "type": []} for m_ in METHODES}
ECART_IAE, ECART_GRANDEURS, ECART_TRAD = 0.0, 0.0, 0.0
NB_TRAD = 0


def rel(a, b):
    """Ecart relatif (absolu si |b| < 1) ; NaN = NaN."""
    if (isinstance(a, float) and math.isnan(a)) or (isinstance(b, float) and math.isnan(b)):
        return 0.0 if (math.isnan(a) and math.isnan(b)) else math.inf
    return abs(a - b) / max(1.0, abs(b))


for m_ in METHODES:
    ns = NS[m_]
    for code in SCENARIOS:
        sc = {s["code"]: s for s in ns["SCENARIOS"] + ns["SCENARIOS_COMPLEMENTAIRES"]}[code]
        reg = FABRIQUES[m_](ns)
        ch = Chrono(reg)
        t0 = time.time()
        sim = ns["simuler"](ch, sc)
        duree = time.time() - t0
        rows = metriques(sim, sc, code)
        d_rows = {(g_, f_, me): v_ for g_, f_, me, v_, _ in rows}
        # Controle 1 : IAE de classement = fichiers (meme expression que les bancs : ecart nul attendu)
        iae = d_rows[("classement", "classement", "IAE")]
        ka_, kb_ = fenetre_classement(code, len(sim["v"]))
        iae_Vs = float(np.sum(np.abs((sim["consigne"] - sim["v"])[ka_:kb_])) * TC)
        if code == "S10":
            att = ATTENDU[m_]["S10_mVs"]
            ecart = abs(iae_Vs * 1e3 - att) / att
        else:
            att = ATTENDU[m_][code] * 1e3
            ecart = abs(iae_Vs - ATTENDU[m_][code]) / ATTENDU[m_][code]
        if abs(iae_Vs * 1e3 - iae) > 1e-12 * iae:
            raise RuntimeError(f"{m_} {code} : IAE incoherente dans le calcul des metriques.")
        ECART_IAE = max(ECART_IAE, ecart)
        if ecart > 1e-12:
            raise RuntimeError(f"{m_} {code} : IAE {iae:.6f} mV.s, fichier {att:.6f} mV.s.")
        # Controle 2 : fonction grandeurs du banc
        o = ns["grandeurs"](sim, sc)
        paires = [(d_rows[("demarrage", "demarrage", "t_montee")], o["t_montee_ms"]),
                  (d_rows[("demarrage", "demarrage", "depassement")], o["depassement_pct"]),
                  (d_rows[("demarrage", "demarrage", "t_etablissement")], o["t_etab_ms"]),
                  (d_rows[("demarrage", "demarrage", "iL_crete")], o["iL_max_dem_A"])]
        for evb in o["evenements"]:
            lab = "ev_%gms" % evb["t_ms"]
            paires += [(d_rows[("evenement", lab, "IAE")], evb["IAE"] * 1e3),
                       (d_rows[("evenement", lab, "e_max")], evb["e_max_V"]),
                       (d_rows[("evenement", lab, "t_retour")],
                        float("nan") if evb["t_retour_ms"] is None else evb["t_retour_ms"]),
                       (d_rows[("evenement", lab, "revenu")], float(evb["revenu"]))]
        if code in ("S2", "S3", "S8a"):
            paires += [(d_rows[("commande", "classement", "dd_moyen")], o["dd_moyen"]),
                       (d_rows[("commande", "classement", "butee")], o["butee_pct"]),
                       (d_rows[("commande", "classement", "dcm")], o["dcm_pct"]),
                       (d_rows[("classement", "classement", "e_max")], o["e_max_V"]),
                       (d_rows[("classement", "classement", "e_eff")], o["e_eff_V"] * 1e3),
                       (d_rows[("classement", "classement", "ISE")], o["ISE"] * 1e3)]
        e_g = max(abs(a - b) / max(abs(b), 1e-300) if not (math.isnan(a) or math.isnan(b))
                  else (0.0 if math.isnan(a) and math.isnan(b) else math.inf) for a, b in paires)
        ECART_GRANDEURS = max(ECART_GRANDEURS, e_g)
        if e_g > 1e-12:
            raise RuntimeError(f"{m_} {code} : metriques differentes de la fonction grandeurs du banc ({e_g:.2e}).")
        # Controle 3 : traduction de Metriques_Simulink.m
        Lt = TR.calculer_metriques(sim["v"], sim["consigne"], sim["d"], sim["iL"], sc["evenements"], code, TE_MATLAB)
        d_t = {(g_, f_, me): v_ for g_, f_, me, v_ in zip(Lt["groupe"], Lt["fenetre"], Lt["metrique"], Lt["valeur"])}
        if set(d_t) != set(d_rows) or len(d_t) != len(Lt["valeur"]):
            raise RuntimeError(f"{m_} {code} : la traduction ne donne pas les memes metriques "
                               f"(en plus : {sorted(set(d_t) - set(d_rows))}, en moins : {sorted(set(d_rows) - set(d_t))}).")
        e_t = max(rel(d_t[k], d_rows[k]) for k in d_rows)
        NB_TRAD += len(d_rows)
        ECART_TRAD = max(ECART_TRAD, e_t)
        if e_t > 1e-9:
            pire = max(d_rows, key=lambda k: rel(d_t[k], d_rows[k]))
            raise RuntimeError(f"{m_} {code} : traduction MATLAB et banc different de {e_t:.2e} sur {pire}.")
        # Temps
        st = stats_temps(ch.t_ns, ch.type, ch.jr is not None)
        TEMPS[m_]["t"] += ch.t_ns
        TEMPS[m_]["type"] += ch.type
        for g_, f_, me, v_, u_ in rows:
            LIGNES.append((code, m_, g_, f_, me, v_, u_))
        for me, (v_, u_) in st.items():
            LIGNES.append((code, m_, "temps_banc", "essai", me, float(v_), u_))
        print(f"  {m_:16s} {code:4s} IAE {iae:8.3f} mV.s (fichier {att:8.3f}, ecart {ecart:.1e}) ; "
              f"grandeurs du banc {e_g:.1e} ; traduction MATLAB {e_t:.1e} ; "
              f"reg.pas moyen {st['appel_moy'][0]:6.2f} us ; simulation {duree:5.1f} s", flush=True)
    st = stats_temps(TEMPS[m_]["t"], TEMPS[m_]["type"], ch.jr is not None)
    for me, (v_, u_) in st.items():
        LIGNES.append(("tous", m_, "temps_banc", "essai", me, float(v_), u_))

for s, m_, g_, f_, me, v_, u_ in LIGNES:
    BASE[TR.cle(s, m_, g_, f_, me)] = v_
print(f"\n  Controle 1 : IAE de classement = fichiers de resultats, ecart relatif maximal {ECART_IAE:.1e}"
      f" ({'nul' if ECART_IAE == 0 else 'NON NUL'}).")
print(f"  Controle 2 : demarrage, evenements, commande = fonction grandeurs du banc, ecart relatif maximal "
      f"{ECART_GRANDEURS:.1e}.")
print(f"  Controle 3 : traduction ligne a ligne de Metriques_Simulink.m, {NB_TRAD} valeurs comparees, ecart "
      f"maximal {ECART_TRAD:.1e} (tolerance 1e-9).")


# %% ETAPE 5 : tableaux, CSV, Markdown

def val(s, m_, g_, f_, me):
    return BASE.get(TR.cle(s, m_, g_, f_, me), float("nan"))


T = ["# Metriques des cinq methodes, banc commun", "",
     "Definitions : en-tete de metriques_banc.py. Grandeurs descriptives ; le classement repose sur l'IAE de la fenetre de classement."]
T = TR.tableaux_communs(T, BASE, SCENARIOS, METHODES)

T += ["", "## Cout de calcul du regulateur seul, banc Python (us par appel de reg.pas ; cinq scenarios reunis)", "",
      "Ordre de grandeur RELATIF entre methodes : Python interprete, un seul processus, pas un materiel embarque. "
      "Le maximum inclut les aleas du systeme.", ""]
T = TR.ligne_md(T, ["Methode", "appel moyen", "appel max", "ordinaire moyen", "fin de fenetre moyen (nombre)",
                    "fin de fenetre max", "avec adaptation moyen (nombre)", "avec adaptation max"])
T = TR.ligne_md(T, ["---"] * 8)
for m_ in METHODES:
    def q(me, fmt="%.2f"):
        return TR.f2(val("tous", m_, "temps_banc", "essai", me), fmt)
    fen = not math.isnan(val("tous", m_, "temps_banc", "essai", "fin_fenetre_moy"))
    T = TR.ligne_md(T, [m_, q("appel_moy"), q("appel_max", "%.1f"),
                        q("ordinaire_moy") if fen else "-",
                        f"{q('fin_fenetre_moy')} ({q('n_fin_fenetre', '%.0f')})" if fen else "-",
                        q("fin_fenetre_max", "%.1f") if fen else "-",
                        f"{q('adaptation_moy')} ({q('n_adaptation', '%.0f')})" if fen else "-",
                        q("adaptation_max", "%.1f") if fen else "-"])
T += ["", "## Cout moyen d'un appel de reg.pas par scenario (us)", ""]
T = TR.ligne_md(T, ["Methode"] + SCENARIOS)
T = TR.ligne_md(T, ["---"] * (len(SCENARIOS) + 1))
for m_ in METHODES:
    T = TR.ligne_md(T, [m_] + [TR.f2(val(s, m_, "temps_banc", "essai", "appel_moy"), "%.2f") for s in SCENARIOS])

with open(os.path.join(DOSSIER, "metriques_banc.csv"), "w", encoding="utf-8", newline="\n") as f:
    f.write("scenario;methode;groupe;fenetre;metrique;valeur;unite\n")
    for s, m_, g_, f_, me, v_, u_ in LIGNES:
        f.write(f"{s};{m_};{g_};{f_};{me};{'NaN' if math.isnan(v_) else format(v_, '.10g')};{u_}\n")
with open(os.path.join(DOSSIER, "metriques_banc.md"), "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(T) + "\n")
titre("ETAPE 5 : tableaux (identiques a metriques_banc.md)")
print("\n".join(T))


print(f"\n  metriques_banc.csv et metriques_banc.md ecrits. Duree totale {time.time() - T_DEBUT:.0f} s.")
