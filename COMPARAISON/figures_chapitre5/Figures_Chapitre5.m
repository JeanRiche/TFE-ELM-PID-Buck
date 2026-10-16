%% FIGURES_CHAPITRE5.m
%
% OBJECTIF
% --------
% Produire les figures de comparaison du chapitre 5 du memoire (scenarios
% S1, S2, S3, S8a et S10), a partir des resultats des
% modeles Simulink de comparaison (<M>_<code>.slx et ZN_<code>.slx, faits
% par Construire_Modeles_Comparaison.m) :
%   Fig5_<code>_A_ZN_ELM.png               Ziegler-Nichols et ELM-PID ;
%   Fig5_<code>_B_ELM_PSO_Fuzzy_PINN.png   ELM-PID, PSO-PID, Fuzzy-PID,
%                                          PINN-PID ;
%   Fig5_S10_gains.png                     gains P, I, D de l'ELM-PID, du
%                                          Fuzzy-PID et du PINN-PID divises
%                                          par ceux de Ziegler-Nichols
%                                          (S10), avec les gains constants
%                                          du PSO-PID.
% Gains du Fuzzy-PID : enregistrement cmp_K de la sortie du bloc
% "Ordonnanceur Flou Zhao" (ordonnanceur_flou_zhao.m), deja les gains
% P, I, D de la forme parallele du bloc PID commun (P = Kp, D = Kd,
% I = Kp^2 / (alpha Kd), eq. 11 de Zhao) : ils
% sont divises par ceux de ZN comme les autres. Ils changent a chaque pas
% avec e et De : la courbe brute est tracee, sans filtrage. Le minimum et
% le maximum de chaque rapport sur l'essai sont affiches dans la console.
%   Fig5_<code>_A_d_ZN_ELM.png             rapport cyclique d, memes
%   Fig5_<code>_B_d_ELM_PSO_Fuzzy_PINN.png methodes, meme mise en page,
%                                          memes zooms et styles que les
%                                          figures de Vout (reglage
%                                          AVEC_RAPPORT_CYCLIQUE).
% Rapport cyclique d : enregistrement cmp_d, pris par
% Construire_Modeles_Comparaison.m sur le signal qui alimente l'entree du
% bloc "PWM Generator (DC-DC)", donc APRES la saturation physique
% [0.01 ; 0.99] du banc commun (D_MIN, D_MAX ; Simuler_Modeles_Comparaison.m
% refuse un d hors de cet intervalle). Dans chaque panneau, une butee
% (0.01 ou 0.99) est tracee en pointilles gris fins seulement si les
% donnees du panneau l'atteignent (a 1e-9 pres) ; elle entre alors dans la
% legende ("saturation (0.01 ; 0.99)"). Ordonnees : donnees du panneau,
% marge de 8 %, rien n'est coupe. La console donne, par scenario et par
% methode, le temps passe en butee (total, avant 30 ms) et le premier et
% le dernier instant en butee.
% Chaque figure est aussi enregistree en .fig. Aucune metrique n'est
% redefinie : le script recalcule seulement l'IAE de la fenetre de
% classement sur les donnees tracees et la compare a metriques_banc.csv
% (controle que les figures viennent des bons modeles).
%
% D'OU VIENNENT LES DONNEES (choix du chemin le plus leger)
% --------------------------------------------------------
% Simuler_Modeles_Comparaison.m (meme script dans les quatre dossiers de
% methode) simule les modeles <M>_<code>.slx et ENREGISTRE les signaux dans
% le dossier de la methode :
%   resultats_Buck_Commun_<M>_comparaison.mat   (ELM_PID, PINN_PID,
%                                                Fuzzy_PID, PSO_PID)
%   resultats_Buck_Commun_comparaison.mat       (Ziegler-Nichols, modeles
%                                                ZN_<code>, dossier PSO_PID)
% variable "resultats" : une entree par scenario, champs code, t, v,
% consigne, d, iL, K (gains P, I, D a chaque pas, methodes adaptatives),
% modele (nom du modele simule). Ce sont exactement les signaux des
% enregistrements cmp_vout et cmp_K des modeles.
% Ce script LIT donc ces fichiers (reglage SOURCE = 'enregistres') :
% quelques secondes au lieu de 15 a 25 minutes de simulation. Si un
% scenario manque dans un fichier, il s'arrete et dit quoi lancer ; avec
% SIMULER_SI_ABSENT = true, il simule alors lui-meme le modele manquant
% (meme lecture que Simuler_Modeles_Comparaison.m : dossier de la methode,
% charger_scenario, sim, enregistrements cmp_vout et cmp_K ; aucun .slx
% n'est enregistre, aucun fichier de resultats n'est ecrit). SOURCE =
% 'simuler' force la simulation de tous les modeles.
% Fichier plus ancien que le modele .slx : avertissement (resultats
% peut-etre perimes, relancer Simuler_Modeles_Comparaison.m).
%
% REGLE DES FENETRES DE ZOOM (la meme pour toutes les methodes et tous les
% resultats)
% ----------------------------------------------------------------------
% Panneau du haut : toute la duree de l'essai (0 a 200 ms).
% Panneaux de zoom, au-dessous :
%   - S1 (aucun evenement) : demarrage, [0 ; 10] ms ;
%   - autres scenarios : pour chaque evenement te du scenario, la fenetre
%     [te - 1 ms ; te + 9 ms]. Les evenements sont lus dans
%     scenario_<code>.mat (variable evenements, ecrite par
%     scenarios_communs.py) et compares a la liste de reference ci-dessous
%     (scenarios_communs.json) :
%       S2, S3 : 50 et 70 ms (debut et fin de la perturbation) ;
%       S8a    : 50, 70, 100, 120 ms ;
%       S10    : 50, 100, 115, 130, 145, 160, 175 ms.
%   - au plus NB_ZOOMS_MAX = 4 panneaux de zoom (grille de deux colonnes) :
%     au-dela (S10 seulement), les quatre PREMIERS evenements sont zoomes
%     (50, 100, 115, 130 ms) ; les suivants (145, 160, 175 ms) ne sont
%     visibles que dans le panneau du haut, ce que la legende du memoire
%     doit dire. Quatre et non trois : les quatre evenements de S8a tiennent
%     tous dans la grille 2 x 2. Les fenetres de zoom sont grisees dans le
%     panneau du haut.
% Ordonnees de chaque panneau : du minimum au maximum des donnees tracees
% DANS ce panneau (toutes les courbes de la figure et la consigne, sur la
% fenetre du panneau), elargi de 8 % de l'etendue de chaque cote, comme
% fixer_ordonnees de Simuler_Modeles_Comparaison.m : aucune courbe ne
% touche le cadre.
%
% STYLES (les memes dans toutes les figures : traits pleins seulement,
% epaisseurs etagees, sans marqueur)
% ------------------------------------------------------------------------
%   Ziegler-Nichols : gris clair [0.65 0.65 0.65], 2.6 pt ;
%   Fuzzy-PID       : vert clair [0.55 0.85 0.45], 3.0 pt ;
%   PSO-PID         : bleu clair [0.35 0.60 0.95], 2.2 pt ;
%   PINN-PID        : noir [0 0 0], 1.5 pt ;
%   ELM-PID         : rouge [0.80 0.10 0.10], 0.9 pt ;
%   consigne        : gris [0.5 0.5 0.5], pointilles ':', 0.8 pt.
% Toutes les courbes des methodes sont en trait plein. Regle de trace : la
% plus epaisse d'abord, la plus fine en dernier (par-dessus). Deux courbes
% confondues se lisent donc comme un trait fin dans une bande plus large
% d'une autre couleur ; la plus fine n'est jamais cachee. Ordre de trace
% dans la figure B : Fuzzy-PID, PSO-PID, PINN-PID, ELM-PID ; figure A :
% ZN puis ELM-PID ; figure des gains : Fuzzy-PID, ZN (droite y = 1), PSO-PID
% (droites), PINN-PID, ELM-PID. Couleur et epaisseur d'une methode sont
% les memes dans toutes les figures (ELM-PID 0.9 pt partout).
% Couleurs : les bandes larges sont claires (vert, bleu, gris), les traits
% fins sont sombres (noir, rouge) ; les deux methodes adaptatives (ELM-PID
% rouge, PINN-PID noir) sont eloignees en teinte et en clarte. Controle
% deuteranopie (matrice de Machado et al., 2009, severite 1, sur RGB
% lineaire ; ecart de couleur dans l'espace OKLab x 100) : ELM-PID vu
% brun-olive [0.51 0.46 0.04], PINN-PID noir, PSO-PID bleu, Fuzzy-PID
% jaune pale, ZN gris. Ecarts des paires d'une meme figure : ELM / PINN
% 57, ELM / PSO 28, ELM / Fuzzy 26, ELM / ZN 20, PINN / PSO 68,
% PINN / Fuzzy 82, PSO / Fuzzy 29 (vision normale : 58, 36, 41, 28, 70,
% 83, 29) ; tous au-dessus de 15, l'epaisseur s'y ajoute.
% Legende : une seule ligne au-dessus de la grille, ordre ELM-PID, PSO-PID,
% Fuzzy-PID, PINN-PID, consigne (figure A : ZN, ELM-PID, consigne ;
% figure des gains : ELM-PID, PSO-PID, Fuzzy-PID, PINN-PID, ZN (= 1)),
% echantillons de trait de 30 points (LONGUEUR_ECHANTILLON) : assez pour
% voir l'epaisseur ; plus longs, cinq entrees ne tiennent plus sur une
% ligne de 16 cm.
% Mise en page : 12 cm de large (taille d'insertion dans le memoire, 0,75 x
% 16 cm), hauteur selon le nombre de rangees (proportions d'une
% figure de 16 cm) ; Times New Roman, etiquettes 10 pt, graduations et
% legendes 9 pt ; pas de titre (la legende de la figure est dans le memoire) ;
% panneaux reperes (a), (b), ... pour la legende du memoire.
% PNG 300 dpi (exportgraphics, sinon print -dpng -r300), .fig, et PDF
% vectoriel (exporter_pdf) sous le nom de l'image du memoire, dans
% pdf_memoire/ (table dans nom_memoire).
%
% CONTROLE DES DONNEES (console)
% ------------------------------
% Pour chaque scenario et chaque methode : IAE de la fenetre de
% classement recalculee sur les donnees tracees, avec la definition de
% calculer_metriques (Metriques_Simulink.m) :
%   S1 : de 0 a 30 ms (IAE de demarrage ; S1 n'a aucun evenement apres
%        30 ms) ;
%   S2, S3, S8a : de 30 ms a la fin ;
%   S10 : de 100 ms a la fin ;
% e = consigne - Vout aux instants k Te, IAE = somme |e| Te (mV.s). Elle
% est comparee a la ligne "classement ; IAE" de
% COMPARAISON/metriques/metriques_banc.csv. Tolerance (comme
% Metriques_Simulink.m) : 1 % ; S10 : 3 % pour l'ELM-PID et le PINN-PID.
% Ziegler-Nichols est compare aussi, a titre indicatif (non controle par
% Metriques_Simulink.m). Un ecart hors tolerance sur une
% autre methode signale un dossier perime ou mal regle (DOSSIER_*) : avec
% ARRETER_SI_ECART = true, le script s'arrete avant de tracer.
%
% UTILISATION
% -----------
% Ce script est dans COMPARAISON/figures_chapitre5 (deux niveaux sous
% la racine des dossiers de methode). Regler DOSSIER_* puis
% lancer. Prerequis : Simuler_Modeles_Comparaison.m lance dans chaque
% dossier de methode (et avec AVEC_ZIEGLER_NICHOLS = true dans PSO_PID),
% scenario_<code>.mat dans les dossiers, metriques_banc.csv dans
% COMPARAISON/metriques. Lecture seule de tous ces fichiers.
% Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

clear; clc;

% --- Reglages ---
% Dossiers de methode, relatifs a RACINE (vide : deux niveaux au-dessus de
% ce script).
RACINE        = '';
DOSSIER_ELM   = 'ELM_PID';
DOSSIER_PINN  = 'PINN_PID';
DOSSIER_PSO   = 'PSO_PID';
DOSSIER_FUZZY = 'FUZZY_PID';
DOSSIER_ZN    = DOSSIER_PSO;                 % modeles ZN_<code> : dossier PSO_PID
SCENARIOS     = {'S1', 'S2', 'S3', 'S8a', 'S10'};
DOSSIER_SORTIE = '';                         % vide : dossier de ce script (figures_chapitre5)
SOURCE = 'enregistres';                      % 'enregistres' : fichiers resultats_*_comparaison.mat ;
                                             % 'simuler' : sim() de chaque modele (15 a 25 min)
SIMULER_SI_ABSENT = false;                   % true : simuler un modele dont le resultat manque
ARRETER_SI_ECART = true;                     % arret si l'IAE s'ecarte du banc hors tolerance
FIGURE_GAINS_AVEC_ZN = true;                 % droite y = 1 (Ziegler-Nichols) dans Fig5_S10_gains
AVEC_RAPPORT_CYCLIQUE = true;                % figures A et B du rapport cyclique d (Fig5_<code>_A_d_..., _B_d_...)
D_BUTEES = [0.01, 0.99];                     % saturation physique devant le PWM (banc commun, D_MIN / D_MAX)
FIGURES_VISIBLES = 'on';                     % 'off' : figures non affichees (enregistrees quand meme)

% Regle des zooms (voir l'en-tete ; ne pas changer d'une figure a l'autre)
ZOOM_S1_MS    = [0, 10];                     % demarrage
ZOOM_AVANT_MS = 1;                           % fenetre [te - 1 ; te + 9] ms
ZOOM_APRES_MS = 9;
NB_ZOOMS_MAX  = 4;                           % au plus 4 zooms (grille 2 x 2)
MARGE_ORDONNEES = 0.08;                      % 8 % de l'etendue de chaque cote
LONGUEUR_ECHANTILLON = 16;                   % longueur des traits de la legende (points)

% Evenements de reference (scenarios_communs.json), en s
EVENEMENTS_REF = struct('S1', zeros(1, 0), 'S2', [0.05, 0.07], 'S3', [0.05, 0.07], ...
                        'S8a', [0.05, 0.07, 0.10, 0.12], ...
                        'S10', [0.05, 0.10, 0.115, 0.13, 0.145, 0.16, 0.175]);

% Mise en forme
LARGEUR_CM = 12;                             % taille d'insertion dans le memoire : 0,75 x 16 cm
LARGEUR_REF_CM = 16;                         % largeur d'origine : les hauteurs gardent ses proportions
POLICE = 'Times New Roman';                  % repli automatique : Helvetica (voir exporter_pdf)
TAILLE_ETIQUETTES = 10;                      % etiquettes d'axes
TAILLE_GRADUATIONS = 9;                      % graduations, legendes, reperes (a), (b), ...
DOSSIER_PDF = 'pdf_memoire';                 % PDF vectoriels, sous les noms du memoire (sous DOSSIER_SORTIE)
RESOLUTION_DPI = 300;
GRIS_CONSIGNE = [0.5 0.5 0.5];
GRIS_FENETRE = [0.93 0.93 0.93];

TE = 1/(22000*10);                           % periode du regulateur et des enregistrements (s)
K_ZN = [0.093910, 301.089, 7.3227e-06];      % P, I, D de Ziegler-Nichols
SEUIL_CONTROLE = 0.01;

% Methodes : nom, dossier, modele source (nom du fichier de resultats),
% prefixe des modeles, couleur, epaisseur (pt, trait plein), tolerance S10.
METH = struct( ...
    'nom',     {'Ziegler-Nichols', 'ELM-PID', 'PSO-PID', 'Fuzzy-PID', 'PINN-PID'}, ...
    'cle',     {'ZN', 'ELM', 'PSO', 'FUZZY', 'PINN'}, ...
    'dossier', {DOSSIER_ZN, DOSSIER_ELM, DOSSIER_PSO, DOSSIER_FUZZY, DOSSIER_PINN}, ...
    'source',  {'Buck_Commun', 'Buck_Commun_ELM_PID', 'Buck_Commun_PSO_PID', 'Buck_Commun_Fuzzy_PID', ...
                'Buck_Commun_PINN_PID'}, ...
    'prefixe', {'ZN', 'ELM_PID', 'PSO_PID', 'Fuzzy_PID', 'PINN_PID'}, ...
    'couleur', {[0.65 0.65 0.65], [0.80 0.10 0.10], [0.35 0.60 0.95], [0.55 0.85 0.45], [0 0 0]}, ...
    'epaisseur', {2.6, 0.9, 2.2, 3.0, 1.5}, ...
    'tol_S10', {0.01, 0.03, 0.01, 0.01, 0.03}, ...
    'controle', {false, true, true, true, true});
FIG_A = {'ZN', 'ELM'};                       % ordre de la legende
FIG_B = {'ELM', 'PSO', 'FUZZY', 'PINN'};

% --- Chemins ---
ICI = fileparts(mfilename('fullpath'));
if isempty(ICI)
    ICI = pwd;
end
if isempty(RACINE)
    RACINE = fullfile(ICI, '..', '..');
end
if isempty(DOSSIER_SORTIE)
    DOSSIER_SORTIE = ICI;
end
fprintf('Racine : %s\nDossiers lus (reglages DOSSIER_* en tete du script) :\n', RACINE);
for m = 1:numel(METH)
    d_ = fullfile(RACINE, METH(m).dossier);
    fprintf('  %-16s : %s\n', METH(m).nom, d_);
    if ~isfolder(d_)
        error('Dossier introuvable : %s (regler DOSSIER_* en tete du script).', d_);
    end
end
fb = fullfile(RACINE, 'COMPARAISON', 'metriques', 'metriques_banc.csv');
if ~isfile(fb)
    error('metriques_banc.csv introuvable (%s) : controle des donnees impossible.', fb);
end
BANC = lire_csv(fb);
fprintf('Source des donnees : %s.\n', SOURCE);

% --- Lecture des donnees et controle ---
D = struct();                                % D.(code).(cle) : t, v, consigne, K
ECARTS = {};
fprintf('\nIAE de la fenetre de classement recalculee sur les donnees tracees (mV.s) :\n');
fprintf('  %-5s %-16s %-14s %10s %10s %8s %6s  %s\n', 'essai', 'methode', 'fenetre', 'figure', 'banc', ...
        'ecart %', 'tol %', 'etat');
for s = 1:numel(SCENARIOS)
    code = SCENARIOS{s};
    for m = 1:numel(METH)
        r = obtenir_essai(RACINE, METH(m), code, SOURCE, SIMULER_SI_ABSENT, TE);
        D.(code).(METH(m).cle) = r;
        [iae, fen] = iae_classement(r.v, r.consigne, code, TE);
        k_ = [code '|' METH(m).nom '|classement|classement|IAE'];
        tol = SEUIL_CONTROLE;
        if strcmp(code, 'S10')
            tol = max(tol, METH(m).tol_S10);
        end
        if isKey(BANC, k_)
            b = BANC(k_);
            ec = iae / b - 1;
            if abs(ec) <= tol
                etat = 'ok';
            elseif METH(m).controle
                etat = 'HORS TOLERANCE';
                ECARTS{end + 1} = sprintf('%s %s : %.3f contre %.3f mV.s au banc (%+.2f %%) - %s', ...
                                          code, METH(m).nom, iae, b, ec * 100, r.fichier); %#ok<SAGROW>
            else
                etat = 'hors tolerance (indicatif)';
            end
        else
            b = NaN;
            ec = NaN;
            etat = 'absent du banc';
        end
        fprintf('  %-5s %-16s %-14s %10.3f %10.3f %+8.2f %6.1f  %s\n', code, METH(m).nom, fen, iae, b, ...
                ec * 100, tol * 100, etat);
    end
end
if ~isempty(ECARTS)
    fprintf('\nATTENTION : donnees qui ne correspondent pas a la methode du banc :\n  %s\n', ...
            strjoin(ECARTS, '\n  '));
    fprintf(['  -> verifier DOSSIER_* (ELM-PID : dossier ELM_PID) et relancer ' ...
             'Simuler_Modeles_Comparaison.m dans le dossier en cause.\n']);
    if ARRETER_SI_ECART
        error('Figures non produites : %d ecart(s) hors tolerance (ARRETER_SI_ECART = true).', numel(ECARTS));
    end
end

% --- Evenements de chaque scenario ---
EV = struct();
for s = 1:numel(SCENARIOS)
    code = SCENARIOS{s};
    EV.(code) = lire_evenements(fullfile(RACINE, DOSSIER_ELM), code, EVENEMENTS_REF);
end

% --- Figures A et B ---
if ~isfolder(DOSSIER_SORTIE)
    mkdir(DOSSIER_SORTIE);
end
FORME = struct('police', POLICE, 'etiq', TAILLE_ETIQUETTES, 'grad', TAILLE_GRADUATIONS, ...
               'largeur', LARGEUR_CM, 'largeur_ref', LARGEUR_REF_CM, 'dossier_pdf', DOSSIER_PDF, 'gris_consigne', GRIS_CONSIGNE, 'gris_fenetre', GRIS_FENETRE, ...
               'marge', MARGE_ORDONNEES, 'dpi', RESOLUTION_DPI, 'visible', FIGURES_VISIBLES, ...
               'echantillon', LONGUEUR_ECHANTILLON);
fprintf('\nFenetres de zoom (regle de l''en-tete) :\n');
for s = 1:numel(SCENARIOS)
    code = SCENARIOS{s};
    [zooms, textes, non_zoomes] = fenetres_zoom(code, EV.(code), ZOOM_S1_MS, ZOOM_AVANT_MS, ZOOM_APRES_MS, ...
                                                NB_ZOOMS_MAX);
    fprintf('  %-4s : %s', code, strjoin(arrayfun(@(i) sprintf('[%g ; %g] ms', zooms(i, 1), zooms(i, 2)), ...
            1:size(zooms, 1), 'UniformOutput', false), ', '));
    if ~isempty(non_zoomes)
        fprintf(' ; evenements non zoomes (panneau du haut seulement) : %s ms', ...
                strjoin(arrayfun(@(x) sprintf('%g', x), non_zoomes, 'UniformOutput', false), ', '));
    end
    fprintf('\n');
    grandeurs_tracees = {'v'};
    if AVEC_RAPPORT_CYCLIQUE
        grandeurs_tracees{end + 1} = 'd';
        resume_butees(D.(code), METH, code, D_BUTEES);
    end
    for gq = grandeurs_tracees
        for f = 1:2
            if strcmp(gq{1}, 'v')
                suffixe = '';
            else
                suffixe = 'd_';
            end
            if f == 1
                cles = FIG_A;
                nom = sprintf('Fig5_%s_A_%sZN_ELM', code, suffixe);
            else
                cles = FIG_B;
                nom = sprintf('Fig5_%s_B_%sELM_PSO_Fuzzy_PINN', code, suffixe);
            end
            fig = figure_vout(D.(code), METH, cles, zooms, textes, FORME, nom, gq{1}, D_BUTEES);
            fichiers = enregistrer(fig, DOSSIER_SORTIE, nom, RESOLUTION_DPI, FORME);
            fprintf('         %s\n', strjoin(fichiers, ' et '));
        end
    end
end

% --- Figure des gains (S10) ---
if any(strcmp(SCENARIOS, 'S10'))
    g_pso = gains_pso(RACINE, DOSSIER_PSO);
    if all(isfinite(g_pso))
        fprintf('\nGains constants du PSO-PID : P = %.5g, I = %.5g, D = %.5g (%.3f / %.3f / %.3f x ZN).\n', ...
                g_pso, g_pso ./ K_ZN);
    end
    fig = figure_gains(D.S10, METH, K_ZN, g_pso, FIGURE_GAINS_AVEC_ZN, FORME);
    fichiers = enregistrer(fig, DOSSIER_SORTIE, 'Fig5_S10_gains', RESOLUTION_DPI, FORME);
    fprintf('Figure des gains : %s\n', strjoin(fichiers, ' et '));
end
fprintf('\nTermine. Figures dans %s.\n', DOSSIER_SORTIE);


%% ===================== Fonctions =====================================

function r = obtenir_essai(racine, M, code, source, simuler_si_absent, Te)
    % Signaux d'un essai : t, v, consigne (colonnes), K (gains a chaque pas,
    % vide pour les gains constants), fichier (origine).
    dossier = fullfile(racine, M.dossier);
    modele = [M.prefixe '_' code];
    f = fullfile(dossier, ['resultats_' M.source '_comparaison.mat']);
    r = [];
    if strcmp(source, 'enregistres')
        if isfile(f)
            S = load(f, 'resultats');
            if isfield(S, 'resultats')
                j = find(strcmp({S.resultats.code}, code), 1);
                if ~isempty(j)
                    e = S.resultats(j);
                    if isfield(e, 'modele') && ~isempty(e.modele) && ~strcmp(e.modele, modele)
                        error('%s : l''entree %s vient du modele %s et non de %s.', f, code, e.modele, modele);
                    end
                    r = normaliser(e, Te);
                    r.fichier = f;
                    fs = fullfile(dossier, [modele '.slx']);
                    if isfile(fs)
                        a = dir(fs);
                        b = dir(f);
                        if a.datenum > b.datenum
                            warning(['%s est plus recent que %s : resultats peut-etre perimes ' ...
                                     '(relancer Simuler_Modeles_Comparaison.m).'], fs, f);
                        end
                    end
                end
            end
        end
        if isempty(r) && ~simuler_si_absent
            error(['%s %s : pas de resultat dans %s.\nLancer Simuler_Modeles_Comparaison.m dans %s ' ...
                   '(MODELES_A_SIMULER = {''%s''}), ou mettre SIMULER_SI_ABSENT = true.'], ...
                  M.nom, code, f, dossier, modele);
        end
    end
    if isempty(r)
        r = simuler_modele(dossier, modele, code, Te);
    end
end

function r = normaliser(e, Te)
    % Colonnes, consigne vectorielle, temps reconstruit s'il manque.
    v = double(e.v(:));
    N = numel(v);
    if isfield(e, 't') && numel(e.t) == N
        t = double(e.t(:));
    else
        t = (0:N - 1)' * Te;
    end
    c = double(e.consigne(:));
    if isscalar(c)
        c = c * ones(N, 1);
    end
    K = [];
    if isfield(e, 'K') && ~isempty(e.K)
        K = double(e.K);
        if size(K, 1) ~= N && size(K, 2) == N
            K = K.';
        end
    end
    d = [];
    if isfield(e, 'd') && numel(e.d) == N
        d = double(e.d(:));
    end
    r = struct('t', t, 'v', v, 'consigne', c, 'K', K, 'd', d, 'fichier', '');
end

function r = simuler_modele(dossier, MDL, code, Te)
    % Meme lecture que Simuler_Modeles_Comparaison.m : dans le dossier de la
    % methode, charger_scenario(code), sim du modele, enregistrements
    % cmp_vout et cmp_K. To Workspace d'origine mis en commentaire EN
    % MEMOIRE ; le modele est ferme sans etre enregistre.
    fs = fullfile(dossier, [MDL '.slx']);
    if ~isfile(fs)
        error('Modele introuvable : %s (lancer Construire_Modeles_Comparaison.m).', fs);
    end
    fprintf('  simulation de %s ...\n', fs);
    ancien = pwd;
    retour = onCleanup(@() cd(ancien));
    cd(dossier);
    sc = charger_scenario(code);
    if bdIsLoaded(MDL)
        if strcmp(get_param(MDL, 'Dirty'), 'on')
            error('%s est ouvert avec des modifications non enregistrees.', MDL);
        end
        close_system(MDL, 0);
    end
    load_system(fs);
    try
        tw = find_system(MDL, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
        noms = get_param(tw, 'VariableName');
        for j = 1:numel(tw)
            if ~startsWith(noms{j}, 'cmp_')
                set_param(tw{j}, 'Commented', 'on');
            end
        end
        sortie = sim(MDL, 'ReturnWorkspaceOutputs', 'on');
    catch ME
        close_system(MDL, 0);
        rethrow(ME);
    end
    close_system(MDL, 0);
    n = numel(sc.t);
    [t, v] = extraire(sortie, 'cmp_vout');
    [~, d] = extraire(sortie, 'cmp_d');
    d = d(1:min(n, end), 1);
    t = t(1:min(n, end));
    v = v(1:min(n, end), 1);
    K = [];
    if any(strcmp(noms, 'cmp_K'))
        [~, K] = extraire(sortie, 'cmp_K');
        K = K(1:min(n, end), :);
    end
    if numel(t) ~= n || max(abs(t - sc.t(:))) > 1e-9
        error('%s : enregistrement non aligne avec le profil.', MDL);
    end
    r = struct('t', t, 'v', v, 'consigne', 100 + double(sc.dvref(:)), 'K', K, 'd', d, 'fichier', [fs ' (sim)']);
    if abs(t(2) - t(1) - Te) > 1e-12
        error('%s : pas d''enregistrement different de Te.', MDL);
    end
end

function [t, M] = extraire(sortie, nom)
    % Enregistrement To Workspace (timeseries ou structure avec temps),
    % une ligne par instant (comme extraire_matrice de
    % Simuler_Modeles_Comparaison.m).
    if ~any(strcmp(sortie.who, nom))
        error('Enregistrement "%s" absent des resultats de simulation.', nom);
    end
    d = sortie.get(nom);
    if isa(d, 'timeseries')
        t = d.Time;
        M = d.Data;
    else
        t = d.time;
        M = d.signals.values;
    end
    t = t(:);
    M = double(squeeze(M));
    if size(M, 1) ~= numel(t)
        M = M.';
    end
end

function [iae, texte] = iae_classement(v, consigne, code, Te)
    % IAE de la fenetre de classement, definition de calculer_metriques
    % (Metriques_Simulink.m) : instants k de [ka ; kb), indice MATLAB k + 1.
    N = numel(v);
    k30 = round(0.03 / Te);
    k100 = round(0.10 / Te);
    if strcmp(code, 'S1')
        ka = 0; kb = k30; texte = '0 a 30 ms';
    elseif strcmp(code, 'S10')
        ka = k100; kb = N; texte = '100 ms a fin';
    else
        ka = k30; kb = N; texte = '30 ms a fin';
    end
    e = consigne(:) - v(:);
    iae = sum(abs(e(ka + 1:kb))) * Te * 1e3;
end

function ev = lire_evenements(dossier, code, ref)
    % Instants des evenements (s) : scenario_<code>.mat, compares a la
    % liste de reference de l'en-tete.
    attendu = ref.(code);
    f = fullfile(dossier, ['scenario_' code '.mat']);
    if isfile(f)
        m = load(f, 'evenements');
        ev = reshape(double(m.evenements), 1, []);
        if numel(ev) ~= numel(attendu) || any(abs(ev - attendu) > 1e-9)
            error('%s : evenements [%s] differents de la reference [%s].', f, num2str(ev), num2str(attendu));
        end
    else
        warning('%s absent : evenements de reference de l''en-tete.', f);
        ev = attendu;
    end
end

function [z, textes, non_zoomes] = fenetres_zoom(code, ev, z_s1, avant, apres, n_max)
    % Regle de l'en-tete. z : une ligne [debut, fin] (ms) par zoom.
    ev_ms = round(ev * 1e3 * 1e6) / 1e6;
    non_zoomes = [];
    if isempty(ev_ms) || strcmp(code, 'S1')
        z = z_s1;
        textes = {['d' lettre(233) 'marrage']};
        return;
    end
    if numel(ev_ms) > n_max
        non_zoomes = ev_ms(n_max + 1:end);
        ev_ms = ev_ms(1:n_max);
    end
    z = [ev_ms(:) - avant, ev_ms(:) + apres];
    textes = arrayfun(@(x) sprintf('%s %g ms', [lettre(233) 'v' lettre(233) 'nement ' lettre(224)], x), ev_ms, ...
                      'UniformOutput', false);
end

function fig = figure_vout(Dc, METH, cles, zooms, textes, F, nom, champ, butees)
    % Panneau du haut (toute la duree) et zooms dessous (grille de une ou
    % deux colonnes). Legende au-dessus, hors des courbes.
    % champ = 'v' : Vout et consigne ; champ = 'd' : rapport cyclique d
    % (entree du PWM) et, dans un panneau seulement si les donnees de ce
    % panneau l'atteignent, la butee de saturation (pointilles gris fins).
    if strcmp(champ, 'd')
        for c = cles
            if isempty(Dc.(c{1}).d)
                error('%s : pas de rapport cyclique d dans les resultats (champ d).', c{1});
            end
        end
    end
    nz = size(zooms, 1);
    nc = 1 + (nz > 1);
    nrz = ceil(nz / nc);
    nr = 1 + nrz;
    hauteur = (6.5 + 4.5 * nrz) * F.largeur / F.largeur_ref;   % cm : 8.25 (1 rangee de zooms), 11.6 (2 rangees)
    fig = figure('Name', nom, 'Color', 'w', 'Units', 'centimeters', 'Visible', F.visible);
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1), max(1, pos(2) - 6), F.largeur, hauteur]);
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [F.largeur, hauteur], ...
        'PaperPosition', [0, 0, F.largeur, hauteur]);
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    idx = cellfun(@(c) find(strcmp({METH.cle}, c)), cles);
    % Ordre de trace : la plus epaisse d'abord, la plus fine par-dessus
    [~, o_] = sort([METH(idx).epaisseur], 'descend');
    ordre_trace = idx(o_);
    lettres = 'abcdefgh';
    ref = Dc.(METH(idx(1)).cle);
    t_ms = ref.t * 1e3;
    fen_tout = [t_ms(1), t_ms(end)];
    fenetres = [fen_tout; zooms];
    ax = gobjects(1, nz + 1);
    h_leg = gobjects(1, numel(idx) + 1);
    butee_vue = false;
    for p = 1:nz + 1
        if p == 1
            ax(p) = nexttile(tl, 1, [1, nc]);
        else
            ligne = 1 + ceil((p - 1) / nc);
            col = mod(p - 2, nc) + 1;
            ax(p) = nexttile(tl, (ligne - 1) * nc + col);
        end
        hold(ax(p), 'on');
        a = fenetres(p, 1);
        b = fenetres(p, 2);
        y_tout = [];
        if p == 1                                  % fenetres de zoom grisees
            for q = 1:nz
                patch(ax(p), zooms(q, [1 2 2 1]), [-1e6 -1e6 1e6 1e6], F.gris_fenetre, ...
                      'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        k = t_ms >= a - 1e-9 & t_ms <= b + 1e-9;
        if strcmp(champ, 'v')
            h_leg(end) = plot(ax(p), t_ms(k), ref.consigne(k), ':', 'Color', F.gris_consigne, 'LineWidth', 0.8);
            y_tout = [y_tout; ref.consigne(k)]; %#ok<AGROW>
        else
            % butees atteintes par les donnees de ce panneau (a 1e-9 pres)
            y_p = [];
            for i = idx
                r = Dc.(METH(i).cle);
                tm = r.t * 1e3;
                y_p = [y_p; r.d(tm >= a - 1e-9 & tm <= b + 1e-9)]; %#ok<AGROW>
            end
            for B = butees
                if any(abs(y_p - B) <= 1e-9)
                    h_leg(end) = plot(ax(p), [a, b], [B, B], ':', 'Color', F.gris_consigne, 'LineWidth', 0.8);
                    butee_vue = true;
                end
            end
        end
        for i = ordre_trace
            r = Dc.(METH(i).cle);
            tm = r.t * 1e3;
            k = tm >= a - 1e-9 & tm <= b + 1e-9;
            h_leg(idx == i) = tracer(ax(p), tm(k), r.(champ)(k), METH(i));
            y_tout = [y_tout; r.(champ)(k)]; %#ok<AGROW>
        end
        xlim(ax(p), [a, b]);
        fixer_ordonnees(ax(p), y_tout, F.marge);
        grid(ax(p), 'on');
        box(ax(p), 'on');
        set(ax(p), 'Layer', 'bottom', 'FontName', F.police, 'FontSize', F.grad);
        if p == 1
            txt = sprintf('(%s) essai complet', lettres(p));
        else
            txt = sprintf('(%s) %s', lettres(p), textes{p - 1});
        end
        title(ax(p), txt, 'FontName', F.police, 'FontSize', F.grad, 'FontWeight', 'normal');
        if p == 1 || p > 1 + (nrz - 1) * nc
            xlabel(ax(p), 'temps (ms)', 'FontName', F.police, 'FontSize', F.etiq);
        end
        if p == 1 || mod(p - 2, nc) == 0
            if strcmp(champ, 'v')
                ylabel(ax(p), 'V_{out} (V)', 'FontName', F.police, 'FontSize', F.etiq);
            else
                ylabel(ax(p), 'rapport cyclique d', 'FontName', F.police, 'FontSize', F.etiq);
            end
        end
    end
    % Objets de legende recrees dans le panneau (a) (donnees NaN, rien de
    % trace) : une legende ne doit citer que des objets de son axe.
    h_leg = gobjects(1, numel(idx) + 1);
    for q = 1:numel(idx)
        h_leg(q) = plot(ax(1), NaN, NaN, '-', 'Color', METH(idx(q)).couleur, 'LineWidth', METH(idx(q)).epaisseur);
    end
    h_leg(end) = plot(ax(1), NaN, NaN, ':', 'Color', F.gris_consigne, 'LineWidth', 0.8);
    if strcmp(champ, 'v')
        noms = [{METH(idx).nom}, {'consigne'}];
    elseif butee_vue
        noms = [{METH(idx).nom}, {strrep(sprintf('butées (%g ; %g)', butees), '.', ',')}];
    else
        noms = {METH(idx).nom};
        delete(h_leg(end));
        h_leg = h_leg(1:end - 1);
    end
    lg = legend(ax(1), h_leg, noms, 'Orientation', 'horizontal', 'Location', 'northoutside', ...
                'FontName', F.police, 'FontSize', F.grad, 'Box', 'off');
    reglages_legende(lg, numel(noms), F);
end

function fig = figure_gains(Dc, METH, K_ZN, g_pso, avec_zn, F)
    % S10 : gains / gains de ZN, un panneau par gain (P, I, D), toute la
    % duree. ELM-PID, Fuzzy-PID et PINN-PID (courbes, enregistrement
    % cmp_K brut, sans filtrage), PSO-PID (droite constante),
    % Ziegler-Nichols (droite y = 1, si demande). Affiche dans la console le
    % minimum et le maximum de chaque rapport sur tout l'essai.
    hauteur = 14 * F.largeur / F.largeur_ref;       % cm : 10.5
    fig = figure('Name', 'Fig5_S10_gains', 'Color', 'w', 'Units', 'centimeters', 'Visible', F.visible);
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1), max(1, pos(2) - 6), F.largeur, hauteur]);
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [F.largeur, hauteur], ...
        'PaperPosition', [0, 0, F.largeur, hauteur]);
    tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    noms_g = {'P', 'I', 'D'};
    adapt = {'ELM', 'FUZZY', 'PINN'};
    pos_leg = [1, 3, 4];                           % ordre de la legende : ELM, PSO, Fuzzy, PINN, ZN
    iadapt = cellfun(@(c) find(strcmp({METH.cle}, c)), adapt);
    ip = find(strcmp({METH.cle}, 'PSO'));
    iz = find(strcmp({METH.cle}, 'ZN'));
    lettres = 'abc';
    h_leg = gobjects(1, 0);
    noms_leg = {};
    for g = 1:3
        ax = nexttile(tl, g);
        hold(ax, 'on');
        y_tout = [];
        t0 = [];
        t1 = [];
        hs = gobjects(1, 5);
        courbes = struct('i', {}, 'pos', {}, 't', {}, 'y', {});
        for q = 1:numel(iadapt)
            i = iadapt(q);
            r = Dc.(METH(i).cle);
            if isempty(r.K)
                warning('%s S10 : pas de gains enregistres (cmp_K) : courbe absente.', METH(i).nom);
                continue;
            end
            tm = r.t * 1e3;
            courbes(end + 1) = struct('i', i, 'pos', pos_leg(q), 't', tm, 'y', r.K(:, g) / K_ZN(g)); %#ok<AGROW>
            t0 = min([t0, tm(1)]);
            t1 = max([t1, tm(end)]);
        end
        if isempty(t0)
            t0 = 0; t1 = 200;
        end
        if all(isfinite(g_pso))
            y = g_pso(g) / K_ZN(g);
            courbes(end + 1) = struct('i', ip, 'pos', 2, 't', [t0; t1], 'y', [y; y]); %#ok<AGROW>
        end
        if avec_zn
            courbes(end + 1) = struct('i', iz, 'pos', 5, 't', [t0; t1], 'y', [1; 1]); %#ok<AGROW>
        end
        % la plus epaisse d'abord, la plus fine par-dessus
        [~, o_] = sort(arrayfun(@(c) METH(c.i).epaisseur, courbes), 'descend');
        for c = courbes(o_)
            hs(c.pos) = tracer(ax, c.t, c.y, METH(c.i));
            y_tout = [y_tout; c.y(:)]; %#ok<AGROW>
        end
        for c = courbes
            fprintf('  S10, %s / %s ZN : %-16s min %8.3f  max %8.3f\n', noms_g{g}, noms_g{g}, METH(c.i).nom, ...
                    min(c.y), max(c.y));
        end
        xlim(ax, [t0, t1]);
        fixer_ordonnees(ax, y_tout, F.marge);
        grid(ax, 'on');
        box(ax, 'on');
        set(ax, 'Layer', 'bottom', 'FontName', F.police, 'FontSize', F.grad);
        ylabel(ax, sprintf('K_%s / K_{%s,ZN}', noms_g{g}, noms_g{g}), 'FontName', F.police, 'FontSize', F.etiq);
        title(ax, sprintf('(%s) gain %s', lettres(g), noms_g{g}), 'FontName', F.police, 'FontSize', F.grad, ...
              'FontWeight', 'normal');
        if g == 3
            xlabel(ax, 'temps (ms)', 'FontName', F.police, 'FontSize', F.etiq);
        end
        if g == 1
            % libelles courts : cinq entrees sur une ligne de 16 cm
            noms4 = {METH(iadapt(1)).nom, METH(ip).nom, METH(iadapt(2)).nom, METH(iadapt(3)).nom, 'ZN (= 1)'};
            ok = arrayfun(@(h) isgraphics(h), hs);
            h_leg = hs(ok);
            noms_leg = noms4(ok);
            ax1 = ax;
        end
    end
    lg = legend(ax1, h_leg, noms_leg, 'Orientation', 'horizontal', 'Location', 'northoutside', ...
                'FontName', F.police, 'FontSize', F.grad, 'Box', 'off');
    reglages_legende(lg, numel(noms_leg), F);
end

function h = tracer(ax, t, y, M)
    % Courbe d'une methode : trait plein, couleur et epaisseur de la
    % methode. Rend l'objet trace (pour la legende).
    h = plot(ax, t(:), y(:), '-', 'Color', M.couleur, 'LineWidth', M.epaisseur);
end

function reglages_legende(lg, n, F)
    % Au-dessus de toute la grille ; une ligne jusqu'a 4 entrees, deux lignes
    % au-dela (sinon la legende deborde des 12 cm de la figure). Traits
    % d'echantillon de F.echantillon points pour montrer l'epaisseur.
    % Proprietes absentes (Octave) : ignorees.
    try
        if n > 4
            lg.NumColumns = ceil(n/2);
        else
            lg.NumColumns = n;
        end
    catch
    end
    try
        lg.ItemTokenSize = [F.echantillon, 18];
    catch
    end
    try
        lg.Layout.Tile = 'north';
    catch
    end
end

function g = gains_pso(racine, dossier)
    % Gains constants du PSO-PID : bloc "PID Controller" du modele
    % PSO_PID_S10.slx si Simulink est la, sinon predictions_banc_pso_pid.json
    % (meme valeurs, ecrites par recherche_pso_pid.py). NaN si rien ne se lit.
    g = NaN(1, 3);
    fs = fullfile(racine, dossier, 'PSO_PID_S10.slx');
    if isfile(fs) && exist('load_system', 'file')
        try
            mdl = 'PSO_PID_S10';
            deja = bdIsLoaded(mdl);
            if ~deja
                load_system(fs);
            end
            noms = {'P', 'I', 'D'};
            for q = 1:3
                g(q) = str2double(get_param([mdl '/PID Controller'], noms{q}));
            end
            if ~deja
                close_system(mdl, 0);
            end
        catch
            g = NaN(1, 3);
        end
    end
    fj = fullfile(racine, dossier, 'predictions_banc_pso_pid.json');
    if isfile(fj)
        j = jsondecode(fileread(fj));
        gj = [j.gains.P, j.gains.I, j.gains.D];
        if all(isfinite(g)) && max(abs(g ./ gj - 1)) > 1e-6
            warning('Gains du PSO-PID : modele [%g %g %g], banc [%g %g %g] : valeurs du modele gardees.', g, gj);
        elseif ~all(isfinite(g))
            g = gj;
        end
    end
    if ~all(isfinite(g))
        warning('Gains du PSO-PID illisibles : pas de droite PSO-PID dans la figure des gains.');
    end
end

function resume_butees(Dc, METH, code, butees)
    % Console : pour chaque methode, temps passe par d sur une butee de
    % saturation (a 1e-9 pres) sur tout l'essai, premier et dernier instant,
    % et la part pendant le demarrage (0 a 30 ms).
    fprintf('         rapport cyclique en butee (%g ou %g), %s :\n', butees, code);
    for m = 1:numel(METH)
        r = Dc.(METH(m).cle);
        if isempty(r.d)
            fprintf('           %-16s : pas de d enregistre\n', METH(m).nom);
            continue;
        end
        te = r.t(2) - r.t(1);
        en_b = abs(r.d - butees(1)) <= 1e-9 | abs(r.d - butees(2)) <= 1e-9;
        dem = r.t < 0.03;
        if ~any(en_b)
            fprintf('           %-16s : jamais\n', METH(m).nom);
        else
            kb = find(en_b);
            fprintf(['           %-16s : %.3f ms au total (%.3f ms avant 30 ms ; %d pas a %g, %d pas a %g), ' ...
                     'de %.3f a %.3f ms\n'], METH(m).nom, sum(en_b) * te * 1e3, sum(en_b & dem) * te * 1e3, ...
                    sum(abs(r.d - butees(1)) <= 1e-9), butees(1), sum(abs(r.d - butees(2)) <= 1e-9), butees(2), ...
                    r.t(kb(1)) * 1e3, r.t(kb(end)) * 1e3);
        end
    end
end

function fixer_ordonnees(ax, y, marge)
    % Copie de Simuler_Modeles_Comparaison.m : [min - m ; max + m], m = marge
    % x etendue ; courbe plate : etendue prise a 10 % de la valeur.
    y = double(y(:));
    y = y(isfinite(y));
    if isempty(y)
        return;
    end
    bas = min(y);
    haut = max(y);
    etendue = haut - bas;
    if etendue <= 0
        etendue = 0.1 * max(abs(haut), 1);
    end
    ylim(ax, [bas - marge * etendue, haut + marge * etendue]);
end

function fichiers = enregistrer(fig, dossier, nom, dpi, F)
    % <nom>.png (exportgraphics, sinon print), <nom>.fig, et le PDF
    % vectoriel sous le nom de l'image du memoire (si le nom y figure).
    fichiers = {fullfile(dossier, [nom '.png']), fullfile(dossier, [nom '.fig'])};
    masquer_barres_outils(fig);
    drawnow;
    try
        exportgraphics(fig, fichiers{1}, 'Resolution', dpi);
    catch
        print(fig, fichiers{1}, '-dpng', sprintf('-r%d', dpi));
    end
    try
        savefig(fig, fichiers{2});
    catch ME
        warning('%s non ecrit (%s).', fichiers{2}, ME.message);
        fichiers = fichiers(1);
    end
    nm = nom_memoire(nom);
    if ~isempty(nm)
        dpdf = fullfile(dossier, F.dossier_pdf);
        if ~isfolder(dpdf)
            mkdir(dpdf);
        end
        fpdf = exporter_pdf(fig, fullfile(dpdf, [nm '.pdf']));
        fichiers{end + 1} = fpdf;
    end
end

function masquer_barres_outils(fig)
    % Barres d'outils des axes masquees avant l'export (sinon MATLAB
    % previent que l'image exportee les montre).
    axes_fig = findall(fig, 'Type', 'axes');
    for j = 1:numel(axes_fig)
        try
            axes_fig(j).Toolbar.Visible = 'off';
        catch
            try
                axtoolbar(axes_fig(j), 'Visible', 'off');
            catch
            end
        end
    end
end

function nm = nom_memoire(nom)
    % Nom (sans extension) de l'image du memoire (MEMOIRE_LATEX/figures/
    % chapitre5) pour un nom du script ; vide si la figure n'est pas dans
    % le memoire (rapport cyclique de S2, S3, S8a : annexe).
    table_noms = {
        'Fig5_S1_A_ZN_ELM',                 'fig5_01_S1_tension_A'
        'Fig5_S1_A_d_ZN_ELM',               'fig5_02_S1_rapport_cyclique_A'
        'Fig5_S2_A_ZN_ELM',                 'fig5_03_S2_tension_A'
        'Fig5_S3_A_ZN_ELM',                 'fig5_04_S3_tension_A'
        'Fig5_S8a_A_ZN_ELM',                'fig5_05_S8a_tension_A'
        'Fig5_S10_A_ZN_ELM',                'fig5_06_S10_tension_A'
        'Fig5_S10_gains',                   'fig5_07_S10_gains'
        'Fig5_S10_A_d_ZN_ELM',              'fig5_08_S10_rapport_cyclique_A'
        'Fig5_S1_B_ELM_PSO_Fuzzy_PINN',     'fig5_09_S1_tension_B'
        'Fig5_S1_B_d_ELM_PSO_Fuzzy_PINN',   'fig5_10_S1_rapport_cyclique_B'
        'Fig5_S2_B_ELM_PSO_Fuzzy_PINN',     'fig5_11_S2_tension_B'
        'Fig5_S3_B_ELM_PSO_Fuzzy_PINN',     'fig5_12_S3_tension_B'
        'Fig5_S8a_B_ELM_PSO_Fuzzy_PINN',    'fig5_13_S8a_tension_B'
        'Fig5_S10_B_ELM_PSO_Fuzzy_PINN',    'fig5_14_S10_tension_B'
        'Fig5_S10_B_d_ELM_PSO_Fuzzy_PINN',  'fig5_15_S10_rapport_cyclique_B'};
    k = find(strcmp(table_noms(:, 1), nom), 1);
    if isempty(k)
        nm = '';
    else
        nm = table_noms{k, 2};
    end
end

function fichier = exporter_pdf(fig, fichier)
    % PDF vectoriel a la taille physique de la figure (12 cm de large :
    % deja fixee dans figure_vout / figure_gains), fond blanc. Ni titre ni
    % reglage de police ici : voir F.police, F.grad, F.etiq. Times New
    % Roman absent de la machine : repli sur Helvetica, avec avertissement.
    fonts = listfonts;
    if ~any(strcmpi(fonts, 'Times New Roman'))
        warning('Times New Roman absente : repli sur Helvetica pour %s.', fichier);
        set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica');
    end
    masquer_barres_outils(fig);
    drawnow;
    if exist('exportgraphics', 'file')
        try
            % Padding 'figure' : garde les 12 cm de large (sinon rognage au contenu).
            exportgraphics(fig, fichier, 'ContentType', 'vector', 'BackgroundColor', 'white', ...
                           'Padding', 'figure');
        catch
            exportgraphics(fig, fichier, 'ContentType', 'vector', 'BackgroundColor', 'white');
        end
    else
        u = get(fig, 'Units');
        set(fig, 'Units', 'centimeters');
        p = get(fig, 'Position');
        set(fig, 'Units', u, 'PaperUnits', 'centimeters', 'PaperSize', p(3:4), ...
            'PaperPosition', [0, 0, p(3:4)], 'InvertHardcopy', 'off');
        print(fig, fichier, '-dpdf', '-vector');
    end
end

function c = lettre(code)
    % Lettre accentuee (code Unicode < 2048) pour le texte des figures.
    % MATLAB : char(code). Octave (caracteres en octets UTF-8) : les deux
    % octets UTF-8 de la lettre.
    if exist('OCTAVE_VERSION', 'builtin')
        c = char([192 + floor(code / 64), 128 + mod(code, 64)]);
    else
        c = char(code);
    end
end

function B = lire_csv(f)
    % metriques_banc.csv -> table cle -> valeur (copie de Metriques_Simulink.m).
    fid = fopen(f, 'r', 'n', 'UTF-8');
    C = textscan(fid, '%s %s %s %s %s %s %s', 'Delimiter', ';', 'Whitespace', '', 'HeaderLines', 1);
    fclose(fid);
    B = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:numel(C{1})
        B([C{1}{i} '|' C{2}{i} '|' C{3}{i} '|' C{4}{i} '|' C{5}{i}]) = str2double(C{6}{i});
    end
end
