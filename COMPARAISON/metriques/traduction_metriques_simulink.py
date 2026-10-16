# =============================================================================
# traduction_metriques_simulink.py
#
# OBJECTIF
#   Traduction Python, ligne a ligne, des fonctions de calcul et de mise en
#   tableaux de Metriques_Simulink.m (calculer_metriques, ajouter,
#   texte_classement, tableaux_communs, fenetres, trier_fenetres,
#   texte_fenetre_rp, f2, valeur, cle). Sans MATLAB, c'est le moyen de
#   controler le script MATLAB : metriques_banc.py applique
#   calculer_metriques ci-dessous aux signaux du banc et verifie qu'elle
#   redonne ses propres metriques (ecrites independamment, en numpy,
#   indices a partir de 0) a 1e-9 pres ; ses tableaux sont produits par
#   tableaux_communs ci-dessous, donc par la meme logique que ceux de
#   MATLAB.
#
# CONVENTION DE TRADUCTION
#   Les indices MATLAB (a partir de 1) sont gardes tels quels ; l'acces
#   MATLAB v(i1:i2) s'ecrit ici m(v, i1, i2) = v[i1-1:i2]. find(x, 1)
#   s'ecrit trouver(x) et renvoie un indice MATLAB (ou None si vide).
#   sum, mean, max, min, abs, diff, sqrt : leurs equivalents numpy (les
#   sommes peuvent differer de MATLAB a 1e-16 pres en relatif, par l'ordre
#   des additions).
#
# UTILISATION
#   Module importe par metriques_banc.py ; il ne s'execute pas seul.
# =============================================================================

import math
import re
import numpy as np


def m(x, i1, i2):
    """x(i1:i2) en MATLAB (indices a partir de 1, i2 compris)."""
    return x[i1 - 1:i2]


def trouver(cond):
    """find(cond, 1) : premier indice MATLAB ou cond est vrai, None si aucun."""
    w = np.flatnonzero(cond)
    return int(w[0]) + 1 if len(w) else None


def trouver_tous(cond):
    """find(cond) : indices MATLAB (a partir de 1)."""
    return np.flatnonzero(cond) + 1


def round_matlab(x):
    """round de MATLAB : demi-entiers arrondis loin de zero."""
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


def g(x):
    """sprintf('%g', x)."""
    return "%g" % x


def ajouter(L, groupe, fenetre, metrique, val, unite):
    L["groupe"].append(groupe)
    L["fenetre"].append(fenetre)
    L["metrique"].append(metrique)
    L["valeur"].append(float(val))
    L["unite"].append(unite)
    return L


def calculer_metriques(v, consigne, d, iL, evenements, code, Te):
    """Traduction de la fonction calculer_metriques de Metriques_Simulink.m."""
    v = np.asarray(v, dtype=float).ravel()
    d = np.asarray(d, dtype=float).ravel()
    iL = np.asarray(iL, dtype=float).ravel()
    N = len(v)
    consigne = np.asarray(consigne, dtype=float).ravel()
    if len(consigne) == 1:
        consigne = consigne * np.ones(N)
    e = consigne - v
    k30 = round_matlab(0.03 / Te)
    k100 = round_matlab(0.10 / Te)
    n5 = round_matlab(0.005 / Te)
    n10 = round_matlab(0.010 / Te)
    evts = [round_matlab(x / Te) for x in np.ravel(evenements)]
    L = {"groupe": [], "fenetre": [], "metrique": [], "valeur": [], "unite": []}
    # 1. Demarrage : [0 ; k_dem)
    k_dem = min([round_matlab(0.05 / Te)] + evts)
    i10 = trouver(v >= 10)
    i90 = trouver(v >= 90)
    if i10 is None or i90 is None:
        t_montee = float("nan")
    else:
        t_montee = (i90 - i10) * Te * 1e3
    L = ajouter(L, "demarrage", "demarrage", "t_montee", t_montee, "ms")
    L = ajouter(L, "demarrage", "demarrage", "depassement", max(0, (np.max(m(v, 1, k_dem)) - 100) / 100 * 100), "%")
    hors = trouver_tous(np.abs(m(e, 1, k_dem)) > 1)
    if len(hors) == 0:
        t_etab = 0
    else:
        t_etab = hors[-1] * Te * 1e3
    L = ajouter(L, "demarrage", "demarrage", "t_etablissement", t_etab, "ms")
    L = ajouter(L, "demarrage", "demarrage", "iL_crete", np.max(m(iL, 1, k_dem)), "A")
    # 2. Fenetre de classement [ka ; kb)
    if code == "S1":
        ka = 0
        kb = k30
    elif code == "S10":
        ka = k100
        kb = N
    else:
        ka = k30
        kb = N
    ew = m(e, ka + 1, kb)                                     # e(idx), idx = (ka + 1):kb
    tau = (np.arange(ka, kb) - ka) * Te
    L = ajouter(L, "classement", "classement", "IAE", np.sum(np.abs(ew)) * Te * 1e3, "mV.s")
    L = ajouter(L, "classement", "classement", "ISE", np.sum(ew ** 2) * Te * 1e3, "V2.ms")
    L = ajouter(L, "classement", "classement", "ITAE", np.sum(tau * np.abs(ew)) * Te * 1e3, "mV.s2")
    L = ajouter(L, "classement", "classement", "e_max", np.max(np.abs(ew)), "V")
    L = ajouter(L, "classement", "classement", "e_eff", np.sqrt(np.mean(ew ** 2)) * 1e3, "mV")
    # 4. Evenements (calcul)
    bornes = evts + [N]
    ne = len(evts)
    ev_iae = np.zeros(ne)
    ev_emax = np.zeros(ne)
    ev_ret = np.zeros(ne)
    ev_rev = np.zeros(ne, dtype=bool)
    for j in range(1, ne + 1):
        a = evts[j - 1]
        b = bornes[j]                                         # bornes(j + 1)
        ef = m(e, a + 1, b)
        ev_iae[j - 1] = np.sum(np.abs(ef)) * Te * 1e3
        ev_emax[j - 1] = np.max(np.abs(ef))
        hors = trouver_tous(np.abs(ef) > 1)
        ev_rev[j - 1] = np.max(np.abs(m(e, max(a, b - n10) + 1, b))) <= 1
        if len(hors) == 0:
            ev_ret[j - 1] = 0
        elif ev_rev[j - 1]:
            ev_ret[j - 1] = hors[-1] * Te * 1e3
        else:
            ev_ret[j - 1] = float("nan")
    n_non = 0
    for j in range(1, ne + 1):
        if evts[j - 1] >= ka and evts[j - 1] < kb and not ev_rev[j - 1]:
            n_non = n_non + 1
    L = ajouter(L, "classement", "classement", "n_non_revenus", n_non, "-")
    # 3. Regime permanent
    prec = [0] + evts
    for j in range(1, ne + 1):
        a = max(prec[j - 1], evts[j - 1] - n5)
        b = evts[j - 1]
        f = "avant_" + g(evts[j - 1] * Te * 1e3) + "ms"
        em = np.mean(m(e, a + 1, b))
        L = ajouter(L, "regime_permanent", f, "e_moy", em * 1e3, "mV")
        L = ajouter(L, "regime_permanent", f, "e_moy_abs", abs(em) * 1e3, "mV")
        L = ajouter(L, "regime_permanent", f, "ondulation_Vout", (np.max(m(v, a + 1, b)) - np.min(m(v, a + 1, b))) * 1e3, "mV")
    a = N - n10
    em = np.mean(m(e, a + 1, N))
    L = ajouter(L, "regime_permanent", "fin", "e_moy", em * 1e3, "mV")
    L = ajouter(L, "regime_permanent", "fin", "e_moy_abs", abs(em) * 1e3, "mV")
    L = ajouter(L, "regime_permanent", "fin", "ondulation_Vout", (np.max(m(v, a + 1, N)) - np.min(m(v, a + 1, N))) * 1e3, "mV")
    # 4. Evenements (suite)
    for j in range(1, ne + 1):
        f = "ev_" + g(evts[j - 1] * Te * 1e3) + "ms"
        L = ajouter(L, "evenement", f, "IAE", ev_iae[j - 1], "mV.s")
        L = ajouter(L, "evenement", f, "e_max", ev_emax[j - 1], "V")
        L = ajouter(L, "evenement", f, "t_retour", ev_ret[j - 1], "ms")
        L = ajouter(L, "evenement", f, "revenu", float(ev_rev[j - 1]), "-")
    # 5. Commande et conduction, fenetre de classement
    dw = m(d, ka + 1, kb)
    L = ajouter(L, "commande", "classement", "dd_moyen", np.mean(np.abs(np.diff(dw))), "-")
    L = ajouter(L, "commande", "classement", "butee", np.mean((dw <= 0.01 + 1e-12) | (dw >= 0.99 - 1e-12)) * 100, "%")
    L = ajouter(L, "commande", "classement", "dcm", np.mean(m(iL, ka + 1, kb) < 0.01) * 100, "%")
    return L


# --- Mise en tableaux (traduction des fonctions du meme nom) ---

def cle(s, m_, g_, f, me):
    return f"{s}|{m_}|{g_}|{f}|{me}"


def valeur(BASE, s, m_, g_, f, me):
    return BASE.get(cle(s, m_, g_, f, me), float("nan"))


def ligne_md(T, cellules):
    T.append("| " + " | ".join(cellules) + " |")
    return T


def f2(x, fmt):
    if isinstance(x, float) and math.isnan(x):
        return "n.d."
    return fmt % x


def texte_classement(noms, vals, etoiles):
    garde = [not math.isnan(x) for x in vals]
    noms = [n for n, k in zip(noms, garde) if k]
    vals = [x for x, k in zip(vals, garde) if k]
    etoiles = [x for x, k in zip(etoiles, garde) if k]
    if not vals:
        return "n.d."
    o = sorted(range(len(vals)), key=lambda i: vals[i])
    vs = [vals[i] for i in o]
    txt = ""
    for i in range(1, len(o) + 1):
        if i > 1:
            if vs[i - 1] / vs[i - 2] - 1 < 0.05:
                txt = txt + " = "
            else:
                txt = txt + " < "
        txt = txt + "%s %.2f%s" % (noms[o[i - 1]], vs[i - 1], etoiles[o[i - 1]])
    return txt


def trier_fenetres(f):
    t = []
    for x in f:
        if x == "fin":
            t.append(math.inf)
        else:
            t.append(float(re.sub(r"^(avant_|ev_)|ms$", "", x)))
    o = sorted(range(len(f)), key=lambda i: t[i])
    return [f[i] for i in o]


def fenetres(BASE, s, METHODES):
    fen_rp, fen_ev = [], []
    cles = sorted(BASE)
    for meth in METHODES:
        pre_rp = f"{s}|{meth}|regime_permanent|"
        pre_ev = f"{s}|{meth}|evenement|"
        for k in cles:
            if k.startswith(pre_rp) and k.endswith("|e_moy"):
                fen_rp.append(k[len(pre_rp):len(k) - len("|e_moy")])
            elif k.startswith(pre_ev) and k.endswith("|IAE"):
                fen_ev.append(k[len(pre_ev):len(k) - len("|IAE")])
        if fen_rp:
            break
    return trier_fenetres(fen_rp), trier_fenetres(fen_ev)


def texte_fenetre_rp(f):
    if f == "fin":
        return "10 dernieres ms"
    return "5 ms " + f.replace("avant_", "avant ").replace("ms", " ms")


FEN = {"S1": "IAE de 0 a 30 ms", "S2": "IAE de 30 ms a la fin", "S3": "IAE de 30 ms a la fin",
       "S8a": "IAE de 30 ms a la fin", "S10": "IAE de 100 ms a la fin"}


def tableaux_communs(T, BASE, SCENARIOS, METHODES):
    # Resume
    T.append("")
    T.append("## Resume : IAE de classement (mV.s) ; erreur moyenne absolue en fin d'essai (mV)")
    T.append("")
    T = ligne_md(T, ["Methode"] + SCENARIOS)
    T = ligne_md(T, ["---"] * (len(SCENARIOS) + 1))
    for meth in METHODES:
        c = [meth]
        for s in SCENARIOS:
            iae = valeur(BASE, s, meth, "classement", "classement", "IAE")
            erp = valeur(BASE, s, meth, "regime_permanent", "fin", "e_moy_abs")
            c.append(f2(iae, "%.2f") + " ; " + f2(erp, "%.2f"))
        T = ligne_md(T, c)
    # Classements
    T.append("")
    T.append('## Classement de chaque scenario (IAE de la fenetre de classement, mV.s ; "=" : moins de 5 % ; '
             '"*" : un evenement de la fenetre non revenu)')
    T.append("")
    for s in SCENARIOS:
        classees = METHODES[1:] if s == "S10" else METHODES
        vals, etoiles = [], []
        for meth in classees:
            vals.append(valeur(BASE, s, meth, "classement", "classement", "IAE"))
            nn = valeur(BASE, s, meth, "classement", "classement", "n_non_revenus")
            etoiles.append("*" if (not math.isnan(nn) and nn > 0) else "")
        T.append("- %s (%s) : %s" % (s, FEN[s], texte_classement(classees, vals, etoiles)))
        if s == "S10":
            zn = valeur(BASE, s, METHODES[0], "classement", "classement", "IAE")
            nn = valeur(BASE, s, METHODES[0], "classement", "classement", "n_non_revenus")
            r = ["%s %s" % (meth, f2(valeur(BASE, s, meth, "classement", "classement", "IAE") / zn, "%.3f"))
                 for meth in METHODES[1:]]
            T.append("  Ziegler-Nichols (reference, hors classement) %s mV.s, %s evenement(s) de la fenetre non "
                     "revenu(s) ; IAE / IAE de ZN : %s" % (f2(zn, "%.2f"), f2(nn, "%.0f"), ", ".join(r)))
    # Par scenario : demarrage, fenetre de classement, commande
    for s in SCENARIOS:
        T.append("")
        T.append("## %s : demarrage, fenetre de classement (%s), commande" % (s, FEN[s]))
        T.append("")
        T = ligne_md(T, ["Methode", "montee 10-90 % (ms)", "depassement (%)", "etabli +-1 V (ms)", "iL crete (A)",
                         "IAE (mV.s)", "ISE (V2.ms)", "ITAE (mV.s2)", "e max (V)", "e eff (mV)", "|dd| moyen",
                         "butee (%)", "DCM (%)"])
        T = ligne_md(T, ["---"] * 13)
        for meth in METHODES:
            T = ligne_md(T, [meth,
                             f2(valeur(BASE, s, meth, "demarrage", "demarrage", "t_montee"), "%.2f"),
                             f2(valeur(BASE, s, meth, "demarrage", "demarrage", "depassement"), "%.2f"),
                             f2(valeur(BASE, s, meth, "demarrage", "demarrage", "t_etablissement"), "%.2f"),
                             f2(valeur(BASE, s, meth, "demarrage", "demarrage", "iL_crete"), "%.2f"),
                             f2(valeur(BASE, s, meth, "classement", "classement", "IAE"), "%.2f"),
                             f2(valeur(BASE, s, meth, "classement", "classement", "ISE"), "%.4g"),
                             f2(valeur(BASE, s, meth, "classement", "classement", "ITAE"), "%.4g"),
                             f2(valeur(BASE, s, meth, "classement", "classement", "e_max"), "%.2f"),
                             f2(valeur(BASE, s, meth, "classement", "classement", "e_eff"), "%.1f"),
                             f2(valeur(BASE, s, meth, "commande", "classement", "dd_moyen"), "%.4f"),
                             f2(valeur(BASE, s, meth, "commande", "classement", "butee"), "%.1f"),
                             f2(valeur(BASE, s, meth, "commande", "classement", "dcm"), "%.1f")])
    # Regime permanent et evenements
    for s in SCENARIOS:
        fen_rp, fen_ev = fenetres(BASE, s, METHODES)
        T.append("")
        T.append("## %s : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete "
                 "de Vout (mV)" % s)
        T.append("")
        T = ligne_md(T, ["Fenetre"] + METHODES)
        T = ligne_md(T, ["---"] * (len(METHODES) + 1))
        for fr in fen_rp:
            c = [texte_fenetre_rp(fr)]
            for meth in METHODES:
                c.append("%s / %s / %s" % (
                    f2(valeur(BASE, s, meth, "regime_permanent", fr, "e_moy"), "%+.2f"),
                    f2(valeur(BASE, s, meth, "regime_permanent", fr, "e_moy_abs"), "%.2f"),
                    f2(valeur(BASE, s, meth, "regime_permanent", fr, "ondulation_Vout"), "%.1f")))
            T = ligne_md(T, c)
        if not fen_ev:
            continue
        T.append("")
        T.append("## %s : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V" % s)
        T.append("")
        T = ligne_md(T, ["Evenement"] + METHODES)
        T = ligne_md(T, ["---"] * (len(METHODES) + 1))
        for fe in fen_ev:
            c = [fe.replace("ev_", "").replace("ms", " ms")]
            for meth in METHODES:
                iae = valeur(BASE, s, meth, "evenement", fe, "IAE")
                if math.isnan(iae):
                    c.append("n.d.")
                    continue
                ret = valeur(BASE, s, meth, "evenement", fe, "t_retour")
                if valeur(BASE, s, meth, "evenement", fe, "revenu") == 0:
                    tr = "pas revenu"
                elif ret == 0:
                    tr = "dans la bande"
                else:
                    tr = "%.2f ms" % ret
                c.append("%.2f / %.2f / %s" % (iae, valeur(BASE, s, meth, "evenement", fe, "e_max"), tr))
            T = ligne_md(T, c)
    return T
