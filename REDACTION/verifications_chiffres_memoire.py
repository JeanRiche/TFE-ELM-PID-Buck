#!/usr/bin/env python3
# =============================================================================
# verifications_chiffres_memoire.py
#
# OBJECTIF
#   Recalculer quatre valeurs citées dans le mémoire qui ne sortent
#   d'aucun autre script du dépôt, à partir des paramètres du circuit et du
#   régulateur que le dépôt fournit, et les comparer au texte :
#     1. le gain statique avec pertes, 196,8 V (soit -1,6 % de 200 V) ;
#     2. la résistance série du condensateur que la somme linéaire des
#        ondulations tolère au pire cas, 0,15 ohm (le calcul sur la forme
#        d'onde, 0,29 ohm, est recalculé en même temps) ;
#     3. le pic de l'action dérivée D.N.dVref sur S2 avec le gain D du PSO-PID,
#        2,45 ;
#     4. les limites de la conduction continue, de 603 à 1 170 ohms.
#
# D'OÙ VIENNENT LES PARAMÈTRES
#   - ELM_PID/banc_commun.py : L, C, Ron, Vf, Rd, fréquence de découpage,
#     bornes du rapport cyclique, gains et coefficient N du PID de
#     Ziegler-Nichols (lus par analyse syntaxique du fichier, sans l'exécuter) ;
#   - PSO_PID/scenarios_communs.json : Vref, tension d'entrée et charge nominales ;
#   - PSO_PID/scenario_S5.mat : tensions d'entrée extrêmes des essais ;
#   - PSO_PID/scenario_S2.mat : saut de consigne de S2 ;
#   - COMPARAISON/banc_pso_pid_resultats.json : gains retenus du PSO-PID.
#   Seules deux grandeurs viennent du cahier des charges du mémoire et d'aucun
#   fichier : la tolérance de 20 % des composants (série E6, section 3.6) et
#   la borne de 0,1 V sur l'ondulation de tension (exigence (iii)).
#
# COMMENT LANCER CE SCRIPT
#   python3 REDACTION/verifications_chiffres_memoire.py
#   Bibliothèques : numpy, scipy. Le code de sortie est 0 si les quatre
#   valeurs coïncident avec le texte, 1 sinon.
# =============================================================================

import ast
import json
import os
import sys

import numpy as np
from scipy.io import loadmat

DEPOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- Grandeurs qui viennent du cahier des charges du mémoire ----------------
TOLERANCE = 0.20          # tolérance relative de L et C (série E6, section 3.6)
ONDULATION_V_MAX = 0.1    # borne de l'ondulation de tension, exigence (iii) (V)


# --- Lecture des paramètres du dépôt -----------------------------------------
def lire_constantes_banc(chemin, noms, ligne_max=160):
    """Lit les constantes du banc Python sans exécuter le fichier.

    Les affectations de tête de fichier (avant ligne_max) sont évaluées dans
    l'ordre, sans fonctions intégrées, avec les constantes déjà lues ; celles
    qui ne se réduisent pas à une expression numérique (chemins, imports)
    sont ignorées.
    """
    with open(chemin, encoding="utf-8") as f:
        arbre = ast.parse(f.read())
    espace = {}
    for noeud in arbre.body:
        if not isinstance(noeud, ast.Assign) or noeud.lineno > ligne_max:
            continue
        cible = noeud.targets[0]
        cibles = cible.elts if isinstance(cible, ast.Tuple) else [cible]
        noms_cibles = [c.id for c in cibles if isinstance(c, ast.Name)]
        valeur = noeud.value
        valeurs = valeur.elts if isinstance(valeur, ast.Tuple) else [valeur]
        if len(noms_cibles) != len(valeurs):
            continue
        for nom, v in zip(noms_cibles, valeurs):
            try:
                expr = ast.Expression(v)
                espace[nom] = eval(compile(expr, chemin, "eval"), {"__builtins__": {}}, dict(espace))
            except Exception:
                continue
    manquants = [n for n in noms if n not in espace]
    if manquants:
        raise RuntimeError(f"constantes introuvables dans {chemin} : {manquants}")
    return {n: espace[n] for n in noms}


def charger_parametres():
    banc = lire_constantes_banc(
        os.path.join(DEPOT, "ELM_PID", "banc_commun.py"),
        ["L_BOB", "C_CONV", "RON", "VF", "RDIODE", "FSW", "D_MIN", "D_MAX",
         "PID_P", "PID_I", "PID_D", "PID_N", "TC"])
    with open(os.path.join(DEPOT, "PSO_PID", "scenarios_communs.json"), encoding="utf-8") as f:
        communs = json.load(f)
    s1 = next(e for e in communs["essais"] if e["code"] == "S1")
    s5 = loadmat(os.path.join(DEPOT, "PSO_PID", "scenario_S5.mat"))
    vin = s5["vin"].ravel()
    s2 = loadmat(os.path.join(DEPOT, "PSO_PID", "scenario_S2.mat"))
    k0 = int(round(float(s2["evenements"].ravel()[0]) / float(s2["Te"].ravel()[0])))
    with open(os.path.join(DEPOT, "COMPARAISON", "banc_pso_pid_resultats.json"), encoding="utf-8") as f:
        pso = json.load(f)
    return {
        "L": banc["L_BOB"], "C": banc["C_CONV"], "Ron": banc["RON"], "Vf": banc["VF"],
        "Rd": banc["RDIODE"], "fs": banc["FSW"], "d_min": banc["D_MIN"], "d_max": banc["D_MAX"],
        "N_zn": banc["PID_N"],
        "Vref": float(communs["Vref"]), "Vin": float(communs["Vin_nominale"]), "R": float(s1["R0"]),
        "Vin_min": float(vin.min()), "Vin_max": float(vin.max()),
        "dVref_S2": float(s2["dvref"].ravel()[k0]),
        "D_pso": float(pso["gains"]["D"]), "N_pso": float(pso["gains"]["N"]),
    }


# --- Calculs -----------------------------------------------------------------
def gain_statique_avec_pertes(p):
    """Gain statique de G_vd(s) avec pertes au point 100 V, 20 A (section 3.9)."""
    Vo = p["Vref"]
    Io = Vo / p["R"]
    d0 = (Vo + p["Vf"] + p["Rd"] * Io) / (p["Vin"] - p["Ron"] * Io + p["Vf"] + p["Rd"] * Io)
    r = d0 * p["Ron"] + (1 - d0) * p["Rd"]
    Ve = p["Vin"] + p["Vf"] - (p["Ron"] - p["Rd"]) * Io
    gain = Ve / (1 + r / p["R"])
    return {"d0": d0, "r": r, "Ve": Ve, "gain": gain, "ecart_pct": 100 * (gain / p["Vin"] - 1)}


def ondulation_totale(p, esr, vin, L, C, n=200001):
    """Ondulation crête à crête de v_C + ESR.i_C sur une période de découpage."""
    Ts = 1 / p["fs"]
    D = p["Vref"] / vin
    di = p["Vref"] * (1 - D) / (p["fs"] * L)
    t = np.linspace(0, Ts, n)
    i = np.where(t < D * Ts, -di / 2 + di * t / (D * Ts), di / 2 - di * (t - D * Ts) / ((1 - D) * Ts))
    vc = np.cumsum(i) * (t[1] - t[0]) / C
    v = vc + esr * i
    return v.max() - v.min()


def limite_esr(p):
    """Résistance série tolérée au pire cas : 240 V, L et C à -20 % (section 3.8.1)."""
    L = (1 - TOLERANCE) * p["L"]
    C = (1 - TOLERANCE) * p["C"]
    vin = p["Vin_max"]
    d_min = p["Vref"] / vin
    di = p["Vref"] * (1 - d_min) / (p["fs"] * L)          # ondulation de courant
    capacitif = 1 / (8 * p["fs"] * C)                       # terme capacitif de (3.19), en ohms
    somme_lineaire = ONDULATION_V_MAX / di - capacitif      # borne de la somme des amplitudes
    # Limite exacte, calculée sur la forme d'onde, par dichotomie.
    bas, haut = 0.0, 1.0
    for _ in range(60):
        milieu = (bas + haut) / 2
        if ondulation_totale(p, milieu, vin, L, C) > ONDULATION_V_MAX:
            haut = milieu
        else:
            bas = milieu
    return {"di": di, "capacitif": capacitif, "somme_lineaire": somme_lineaire,
            "exacte": bas,
            "ond_029": ondulation_totale(p, 0.29, vin, L, C),
            "ond_030": ondulation_totale(p, 0.30, vin, L, C)}


def pic_derivee(p):
    """Pic D.N.dVref de l'incrément dérivé de la loi incrémentale (section 4.7.2)."""
    pic = p["D_pso"] * p["N_pso"] * p["dVref_S2"]
    return {"pic": pic, "plage_d": p["d_max"] - p["d_min"], "N_coherent": abs(p["N_pso"] - p["N_zn"]) < 1e-9}


def limites_conduction(p):
    """R < 2 L fs / (1 - d0) avec d0 = Vo/Vin, sans pertes (section 3.8.1)."""
    def rcrit(vin, L):
        return 2 * L * p["fs"] / (1 - p["Vref"] / vin)
    return {"nominal": rcrit(p["Vin"], p["L"]),
            "vin_min": rcrit(p["Vin_min"], p["L"]),
            "pire_cas": rcrit(p["Vin_max"], (1 - TOLERANCE) * p["L"])}


# --- Affichage et comparaison au texte ----------------------------------------
def ligne(nom, calcul, texte, tolerance, unite, section, decimales):
    ecart = calcul - texte
    ok = abs(ecart) <= tolerance
    fmt = f"{{:.{decimales}f}}"
    print(f"  {nom:<46} calcul {fmt.format(calcul):>10} {unite:<3} | texte {fmt.format(texte):>10} | "
          f"{'coïncide' if ok else 'ÉCART':<8} | {section}")
    return ok


def main():
    p = charger_parametres()
    print("Paramètres lus dans le dépôt : "
          f"L = {p['L']*1e3:g} mH, C = {p['C']*1e6:g} uF, fs = {p['fs']/1e3:g} kHz, Ron = {p['Ron']} ohm, "
          f"Vf = {p['Vf']} V, Rd = {p['Rd']*1e3:g} mohm, Vin = {p['Vin']:g} V "
          f"(essais de {p['Vin_min']:g} à {p['Vin_max']:g} V), R = {p['R']:g} ohm, Vref = {p['Vref']:g} V\n")
    resultats = []

    g = gain_statique_avec_pertes(p)
    print("1. Gain statique avec pertes (point 100 V, 20 A)")
    print(f"     d0 = {g['d0']:.4f}, r = {g['r']:.4f} ohm, Ve = {g['Ve']:.3f} V")
    resultats.append(ligne("gain statique de G_vd(s)", g["gain"], 196.8, 0.05, "V", "section 3.9 (sec:base-modele-nominal)", 3))
    resultats.append(ligne("écart au gain sans pertes de 200 V", g["ecart_pct"], -1.6, 0.05, "%", "section 3.9 (sec:base-modele-nominal)", 2))
    resultats.append(ligne("rapport cyclique d0 pour 100 V, 20 A", g["d0"], 0.5071, 0.00005, "", "section 3.9 (sec:base-modele-nominal)", 4))
    resultats.append(ligne("résistance série moyenne r", g["r"], 0.05, 0.005, "ohm", "section 3.9 (sec:base-modele-nominal)", 4))
    resultats.append(ligne("tension équivalente Ve", g["Ve"], 198.8, 0.05, "V", "section 3.9 (sec:base-modele-nominal)", 3))

    e = limite_esr(p)
    print("\n2. Résistance série du condensateur au pire cas (240 V, L et C à -20 %)")
    print(f"     ondulation de courant = {e['di']:.4f} A, 1/(8 fs C) = {e['capacitif']:.4f} ohm")
    resultats.append(ligne("limite par la somme linéaire (3.19)", e["somme_lineaire"], 0.15, 0.005, "ohm", "section 3.8.1 (sec:base-verification)", 4))
    # Le texte retient 0,29 ohm, valeur arrondie par défaut : à 0,30 ohm l'ondulation dépasse 0,1 V de 0,1 mV.
    print(f"     limite exacte sur la forme d'onde = {e['exacte']:.4f} ohm (le texte l'arrondit par défaut à 0,01 ohm près)")
    resultats.append(ligne("limite sur la forme d'onde, arrondie par défaut", np.floor(e["exacte"] * 100) / 100, 0.29, 0.0005, "ohm", "section 3.8.1 (sec:base-verification)", 4))
    resultats.append(ligne("ondulation totale à 0,29 ohm", e["ond_029"], 0.097, 0.0005, "V", "section 3.8.1 (sec:base-verification)", 4))
    resultats.append(ligne("ondulation totale à 0,30 ohm", e["ond_030"], 0.100, 0.0005, "V", "section 3.8.1 (sec:base-verification)", 4))

    d = pic_derivee(p)
    print("\n3. Pic de l'action dérivée sur S2 avec le D du PSO-PID")
    print(f"     D = {p['D_pso']:.6e}, N = {p['N_pso']:.1f} rad/s (N du bloc commun : "
          f"{'identique' if d['N_coherent'] else 'DIFFÉRENT'}), saut de consigne = {p['dVref_S2']:g} V, "
          f"plage du rapport cyclique = {d['plage_d']:.2f}")
    resultats.append(ligne("pic D.N.dVref", d["pic"], 2.45, 0.005, "", "section 4.7.2", 4))
    resultats.append(d["N_coherent"] and d["pic"] > d["plage_d"])

    c = limites_conduction(p)
    print("\n4. Limites de la conduction continue, R < 2 L fs / (1 - d0)")
    resultats.append(ligne("200 V, L nominale", c["nominal"], 880, 5, "ohm", "section 3.8.1 (sec:base-verification)", 1))
    resultats.append(ligne("160 V, L nominale (arrondi à la dizaine)", c["vin_min"], 1170, 5, "ohm", "section 3.8.1 et 6.2.2 (sec:persp-dynamique-non-lineaire-sans)", 1))
    resultats.append(ligne("240 V, L à -20 % (pire cas)", c["pire_cas"], 603, 0.5, "ohm", "section 3.8.1 et 6.2.2 (sec:persp-dynamique-non-lineaire-sans)", 1))

    ok = all(resultats)
    print("\nBilan : " + ("toutes les valeurs coïncident avec le texte du mémoire." if ok
                        else "au moins une valeur s'écarte du texte du mémoire."))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
