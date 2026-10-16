%% METRIQUES_SIMULINK.m
%
% OBJECTIF
% --------
% Calculer, sur les resultats SIMULINK des cinq methodes, les metriques
% descriptives definies dans l'en-tete de metriques_banc.py, avec
% EXACTEMENT les memes definitions que ce script :
%   1. demarrage : temps de montee 10-90 %, depassement, temps
%      d'etablissement +-1 V, courant crete iL ;
%   2. fenetre de classement (S1 : 0 a 30 ms ; S2, S3, S8a : 30 ms a la
%      fin ; S10 : 100 ms a la fin) : IAE, ISE, ITAE (temps compte depuis
%      le debut de la fenetre), ecart maximal, ecart efficace ;
%   3. erreur en regime permanent : moyenne de (consigne - Vout) sur les
%      5 ms avant chaque evenement et sur les 10 dernieres ms de l'essai
%      (signee et absolue), ondulation crete a crete de Vout ;
%   4. par evenement : IAE, ecart maximal, temps de retour dans +-1 V ;
%   5. activite de la commande : |Delta d| moyen, temps en butee,
%      conduction discontinue (fenetre de classement) ;
%   6b. temps de calcul Simulink par periode du regulateur (duree_calcul_s
%      et us_par_pas des fichiers de resultats). Ce n'est PAS le cout du
%      regulateur : l'execution interpretee des blocs MATLAB System et le
%      circuit de puissance au pas du powergui dominent.
% Les metriques ne changent aucun classement : le classement affiche est
% celui de l'IAE de la fenetre de classement (regle des 5 %), S10 avec
% quatre methodes classees et Ziegler-Nichols en reference.
%
% FICHIERS LUS (aucun n'est modifie)
% ---------------------------------
% Ce script est dans COMPARAISON/metriques. Il lit, a partir de la racine
% du depot (noms des dossiers : reglages DOSSIER_* en tete du script ; ils
% sont affiches au lancement) :
%   ../../PSO_PID/resultats_Buck_Commun_PSO_PID.mat      et ..._S10.mat
%   ../../FUZZY_PID/resultats_Buck_Commun_Fuzzy_PID.mat  et ..._S10.mat
%   ../../ELM_PID/resultats_Buck_Commun_ELM_PID.mat      et ..._S10.mat
%   ../../PINN_PID/resultats_Buck_Commun_PINN_PID.mat    et ..._S10.mat
%   Ziegler-Nichols : resultats_Buck_Commun.mat et resultats_Buck_Commun_S10.mat,
%   cherches dans les quatre dossiers (ecrits par Simuler_ZN.m, dossier
%   PSO_PID). Le premier trouve est pris.
%   Pour S10, le fichier ..._S10.mat (Simuler_<M>.m avec ESSAIS_A_SIMULER =
%   {'S10'}) passe avant le fichier principal s'il contient aussi S10.
%   Les instants des evenements sont lus dans scenario_<code>.mat du meme
%   dossier.
%   Facultatif : metriques_banc.csv (ce dossier, ecrit par metriques_banc.py)
%   pour comparer l'IAE de classement Simulink a celle du banc.
%   Avec le reglage RESULTATS = 'comparaison', le script lit a la place les
%   fichiers resultats_<modele>_comparaison.mat ecrits par
%   Simuler_Modeles_Comparaison.m (modeles lances d'un clic,
%   <PREFIXE>_<code>.slx, de Construire_Modeles_Comparaison.m), un fichier
%   par methode pour les cinq scenarios, au meme format ; les sorties sont
%   alors metriques_simulink_comparaison.csv et .md. Les definitions et le
%   calcul ne changent pas.
% Un fichier absent est signale ; le script continue avec ce qu'il trouve.
%
% CONTROLE DE CORRESPONDANCE
% --------------------------
% Les fichiers de resultats ne contiennent pas d'identifiant du regulateur
% simule. Le controle compare donc, pour chaque fichier, l'IAE de
% classement de chaque scenario a celle du banc (metriques_banc.csv). Si un
% scenario s'en ecarte de plus de SEUIL_CONTROLE (1 %), le fichier entier
% est REFUSE (aucune de ses valeurs n'est utilisee) avec le message "ce
% fichier ne correspond pas a la methode du banc". Exceptions : S10,
% tolerance plus large si elle est documentee (3 % ELM-PID et PINN-PID) ;
% S7a (quantification, ecarts de 1 a 6 %) n'est pas controle.
% Ziegler-Nichols n'est pas controle. Simulink retrouve le banc a 0.13 %
% pres (Simuler_Modeles_Comparaison.m, S1 a S10) : le seuil de 1 % ne
% refuse pas ces fichiers. metriques_banc.csv est donc obligatoire.
%
% CE QUE PRODUIT CE SCRIPT (dans ce dossier)
% ------------------------------------------
%   - dans la console : les dossiers lus, les fichiers refuses par le
%     controle, les memes tableaux que metriques_banc.py (format Markdown)
%     et, si metriques_banc.csv est la, l'ecart Simulink - banc sur l'IAE
%     de classement (S10 : tolerance 1 % PSO-PID et Fuzzy-PID, 3 % ELM-PID
%     et PINN-PID) ;
%   - metriques_simulink.csv : une ligne par metrique (scenario ; methode ;
%     groupe ; fenetre ; metrique ; valeur ; unite), separateur ';' ;
%   - metriques_simulink.md : les tableaux, prets a copier.
%   (RESULTATS = 'comparaison' : metriques_simulink_comparaison.csv et .md)
%
% CONTROLE SANS MATLAB
% --------------------
% La fonction calculer_metriques (en bas du script) est traduite ligne a
% ligne dans traduction_metriques_simulink.py. metriques_banc.py applique
% cette traduction aux signaux du banc et verifie qu'elle redonne ses
% propres metriques a 1e-9 pres.
%
% Duree : quelques secondes (lecture des fichiers .mat). Compatible R2024a.
% ---------------------------------------------------------------------

clear; clc;

% --- Reglages ---
RESULTATS = 'Simuler';      % 'Simuler' : fichiers de Simuler_<M>.m (resultats_<modele>.mat et _S10.mat) ;
                            % 'comparaison' : fichiers de Simuler_Modeles_Comparaison.m
                            % (resultats_<modele>_comparaison.mat)
% Dossiers de methode, relatifs a la racine du depot (deux niveaux au-dessus
% de ce dossier).
DOSSIER_PSO   = 'PSO_PID';
DOSSIER_FUZZY = 'FUZZY_PID';
DOSSIER_ELM   = 'ELM_PID';
DOSSIER_PINN  = 'PINN_PID';
% Controle de correspondance (voir l'en-tete) : ecart relatif maximal entre l'IAE
% de classement Simulink et celle du banc (metriques_banc.csv).
SEUIL_CONTROLE = 0.01;
SCENARIOS_SANS_CONTROLE = {'S7a'};                         % quantification : ecarts documentes de 1 a 6 %
CONTROLE_METHODE = [false, true, true, true, true];        % Ziegler-Nichols : non controle
SCENARIOS = {'S1', 'S2', 'S3', 'S8a', 'S10'};
METHODES  = {'Ziegler-Nichols', 'PSO-PID', 'Fuzzy-PID', 'ELM-PID', 'PINN-PID'};
DOSSIERS  = {{DOSSIER_PSO, DOSSIER_FUZZY, DOSSIER_ELM, DOSSIER_PINN}, {DOSSIER_PSO}, {DOSSIER_FUZZY}, ...
             {DOSSIER_ELM}, {DOSSIER_PINN}};               % ou chercher les resultats
MODELES   = {'Buck_Commun', 'Buck_Commun_PSO_PID', 'Buck_Commun_Fuzzy_PID', ...
             'Buck_Commun_ELM_PID', 'Buck_Commun_PINN_PID'};
TOLERANCE_S10 = [NaN, 0.01, 0.01, 0.03, 0.03];             % IAE de 100 ms a la fin, Simulink / banc
TE = 1/(22000*10);                                         % periode du regulateur et des enregistrements (s)

ICI = fileparts(mfilename('fullpath'));
if isempty(ICI)
    ICI = pwd;
end
RACINE = fullfile(ICI, '..', '..');
switch RESULTATS
    case 'Simuler'
        SUFFIXES = {'', '_S10'};                           % le fichier _S10 passe avant pour S10
        SORTIE = 'metriques_simulink';
    case 'comparaison'
        SUFFIXES = {'_comparaison'};
        SORTIE = 'metriques_simulink_comparaison';
    otherwise
        error('RESULTATS = ''%s'' : choisir ''Simuler'' ou ''comparaison''.', RESULTATS);
end
fprintf('Resultats lus : %s (fichiers resultats_<modele>%s.mat).\n', RESULTATS, strjoin(SUFFIXES, '.mat ou '));
fprintf('Dossiers lus (reglages DOSSIER_* en tete du script) :\n');
for im = 1:numel(METHODES)
    chemins = cellfun(@(x) fullfile(RACINE, x), DOSSIERS{im}, 'UniformOutput', false);
    fprintf('  %-16s : %s\n', METHODES{im}, strjoin(chemins, ' ; '));
    for q = 1:numel(chemins)
        if ~isfolder(chemins{q})
            fprintf('  ATTENTION : dossier introuvable : %s (regler DOSSIER_* en tete du script)\n', chemins{q});
        end
    end
end
fb = fullfile(ICI, 'metriques_banc.csv');
if ~isfile(fb)
    error(['metriques_banc.csv absent de %s : le controle de correspondance (fichier de la methode retenue ' ...
           'de chaque methode) est impossible. Reprendre le fichier du depot.'], ICI);
end
BANC = lire_csv(fb);

% --- Lecture des resultats et calcul des metriques ---
LIGNES = cell(0, 7);                                       % scenario, methode, groupe, fenetre, metrique, valeur, unite
ABSENTS = {};
REFUSES = {};                                              % fichiers refuses par le controle de correspondance
for im = 1:numel(METHODES)
    E = lire_resultats(RACINE, DOSSIERS{im}, MODELES{im}, SUFFIXES);
    if CONTROLE_METHODE(im)
        [E, refus] = controler_correspondance(E, METHODES{im}, BANC, SEUIL_CONTROLE, TOLERANCE_S10(im), ...
                                       SCENARIOS_SANS_CONTROLE, TE);
        REFUSES = [REFUSES, refus]; %#ok<AGROW>
    end
    if isempty(E)
        fprintf('%-16s : aucun fichier resultats_%s*.mat dans %s.\n', METHODES{im}, MODELES{im}, ...
                strjoin(DOSSIERS{im}, ', '));
    end
    for is = 1:numel(SCENARIOS)
        code = SCENARIOS{is};
        en = choisir(E, code);
        if isempty(en)
            ABSENTS{end + 1} = sprintf('%s %s', METHODES{im}, code); %#ok<SAGROW>
            continue;
        end
        r = en.r;
        evts = lire_evenements(en);
        if ~isequal(numel(r.v), numel(r.d), numel(r.iL))
            fprintf('  ATTENTION : %s %s, signaux de longueurs differentes dans %s : ignore.\n', ...
                    METHODES{im}, code, en.fichier);
            ABSENTS{end + 1} = sprintf('%s %s (signaux incoherents)', METHODES{im}, code); %#ok<SAGROW>
            continue;
        end
        L = calculer_metriques(r.v, r.consigne, r.d, r.iL, evts, code, TE);
        for q = 1:numel(L.valeur)
            LIGNES(end + 1, :) = {code, METHODES{im}, L.groupe{q}, L.fenetre{q}, L.metrique{q}, ...
                                  L.valeur(q), L.unite{q}}; %#ok<SAGROW>
        end
        LIGNES(end + 1, :) = {code, METHODES{im}, 'temps_simulink', 'essai', 'duree_calcul', ...
                              double(r.duree_calcul_s), 's'}; %#ok<SAGROW>
        LIGNES(end + 1, :) = {code, METHODES{im}, 'temps_simulink', 'essai', 'par_periode_Tc', ...
                              double(r.us_par_pas), 'us'}; %#ok<SAGROW>
        fprintf('%-16s %-4s : lu dans %s\n', METHODES{im}, code, en.fichier);
    end
end
if ~isempty(REFUSES)
    fprintf(['\nFICHIERS REFUSES (ne correspondent pas a la methode du banc ; leurs scenarios sont absents) :' ...
             '\n  %s\n'], strjoin(REFUSES, [newline '  ']));
end
if ~isempty(ABSENTS)
    fprintf('\nResultats Simulink ABSENTS (non calcules) :\n  %s\n', strjoin(ABSENTS, [newline '  ']));
    if strcmp(RESULTATS, 'Simuler') && any(contains(ABSENTS, 'S10'))
        fprintf(['  S10 : si S10 a ete simule avec les modeles de comparaison (<PREFIXE>_S10.slx, ' ...
                 'Simuler_Modeles_Comparaison.m), relancer avec RESULTATS = ''comparaison''.\n']);
    end
end
BASE = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:size(LIGNES, 1)
    BASE(cle(LIGNES{i, 1}, LIGNES{i, 2}, LIGNES{i, 3}, LIGNES{i, 4}, LIGNES{i, 5})) = LIGNES{i, 6};
end

% --- CSV ---
fcsv = fullfile(ICI, [SORTIE '.csv']);
fid = fopen(fcsv, 'w', 'n', 'UTF-8');
fprintf(fid, 'scenario;methode;groupe;fenetre;metrique;valeur;unite\n');
for i = 1:size(LIGNES, 1)
    fprintf(fid, '%s;%s;%s;%s;%s;%s;%s\n', LIGNES{i, 1:5}, nombre_csv(LIGNES{i, 6}), LIGNES{i, 7});
end
fclose(fid);

% --- Tableaux (Markdown, affiches aussi dans la console) ---
T = {};
if strcmp(RESULTATS, 'comparaison')
    T{end + 1} = '# Metriques des cinq methodes, resultats Simulink des modeles de comparaison (<PREFIXE>_<code>.slx)';
else
    T{end + 1} = '# Metriques des cinq methodes, resultats Simulink';
end
T{end + 1} = '';
T{end + 1} = 'Definitions : definitions_metriques.txt. Grandeurs descriptives ; le classement reste celui des IAE deja fixees.';
if ~isempty(REFUSES)
    T{end + 1} = '';
    T{end + 1} = ['Fichiers refuses par le controle de correspondance (ne correspondent pas a la methode du banc) : ' ...
                  strjoin(REFUSES, ' ; ')];
end
if ~isempty(ABSENTS)
    T{end + 1} = '';
    T{end + 1} = ['Resultats absents (n.d.) : ' strjoin(ABSENTS, ', ') '.'];
end
T = tableaux_communs(T, BASE, SCENARIOS, METHODES);

% Temps de calcul Simulink (6b)
T{end + 1} = '';
T{end + 1} = '## Temps de calcul Simulink par periode du regulateur (us ; duree de sim() entre parentheses)';
T{end + 1} = '';
T{end + 1} = 'Lu dans les resultats (us_par_pas = duree de sim() / nombre de periodes Tc). Ce n''est pas le cout du regulateur : l''execution interpretee et le circuit de puissance au pas du powergui dominent.';
T{end + 1} = '';
T = ligne_md(T, [{'Methode'}, SCENARIOS]);
T = ligne_md(T, repmat({'---'}, 1, numel(SCENARIOS) + 1));
for im = 1:numel(METHODES)
    c = {METHODES{im}};
    for is = 1:numel(SCENARIOS)
        us = valeur(BASE, SCENARIOS{is}, METHODES{im}, 'temps_simulink', 'essai', 'par_periode_Tc');
        du = valeur(BASE, SCENARIOS{is}, METHODES{im}, 'temps_simulink', 'essai', 'duree_calcul');
        if isnan(us)
            c{end + 1} = 'n.d.'; %#ok<SAGROW>
        else
            c{end + 1} = sprintf('%.1f (%.0f s)', us, du); %#ok<SAGROW>
        end
    end
    T = ligne_md(T, c);
end

% Ecart Simulink - banc sur l'IAE de classement (BANC lu au debut)
T{end + 1} = '';
T{end + 1} = '## IAE de classement : Simulink / banc (ecart)';
T{end + 1} = '';
if ~isfile(fb)
    T{end + 1} = 'metriques_banc.csv absent de ce dossier : pas de comparaison au banc.';
else
    T{end + 1} = 'S10 : tolerance de 1 % pour PSO-PID et Fuzzy-PID, 3 % pour ELM-PID et PINN-PID ; pas de tolerance pour Ziegler-Nichols.';
    T{end + 1} = '';
    T = ligne_md(T, [{'Methode'}, SCENARIOS]);
    T = ligne_md(T, repmat({'---'}, 1, numel(SCENARIOS) + 1));
    for im = 1:numel(METHODES)
        c = {METHODES{im}};
        for is = 1:numel(SCENARIOS)
            k = cle(SCENARIOS{is}, METHODES{im}, 'classement', 'classement', 'IAE');
            s = valeur(BASE, SCENARIOS{is}, METHODES{im}, 'classement', 'classement', 'IAE');
            if isKey(BANC, k) && ~isnan(s)
                b = BANC(k);
                txt = sprintf('%.2f / %.2f (%+.2f %%)', s, b, (s / b - 1) * 100);
                if strcmp(SCENARIOS{is}, 'S10') && ~isnan(TOLERANCE_S10(im))
                    if abs(s / b - 1) <= TOLERANCE_S10(im)
                        txt = [txt ' conforme']; %#ok<AGROW>
                    else
                        txt = [txt ' HORS TOLERANCE']; %#ok<AGROW>
                    end
                end
                c{end + 1} = txt; %#ok<SAGROW>
            else
                c{end + 1} = 'n.d.'; %#ok<SAGROW>
            end
        end
        T = ligne_md(T, c);
    end
end

fmd = fullfile(ICI, [SORTIE '.md']);
fid = fopen(fmd, 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', T{:});
fclose(fid);
fprintf('\n');
fprintf('%s\n', T{:});
fprintf('\n%s et %s ecrits.\n', fcsv, fmd);


%% ===================== Fonctions =====================================

function L = calculer_metriques(v, consigne, d, iL, evenements, code, Te)
    % Metriques de definitions_metriques.txt pour un essai. Les numeros
    % d'instants k commencent a 0 (comme en Python) ; l'indice MATLAB de
    % l'instant k est k + 1. Traduite ligne a ligne dans
    % traduction_metriques_simulink.py.
    v = double(v(:));
    d = double(d(:));
    iL = double(iL(:));
    N = numel(v);
    consigne = double(consigne(:));
    if numel(consigne) == 1
        consigne = consigne * ones(N, 1);
    end
    e = consigne - v;
    k30 = round(0.03 / Te);
    k100 = round(0.10 / Te);
    n5 = round(0.005 / Te);                              % 5 ms : 1 100 instants
    n10 = round(0.010 / Te);                             % 10 ms : 2 200 instants
    evts = round(reshape(evenements, 1, []) / Te);       % instants des evenements (k)
    L.groupe = {};
    L.fenetre = {};
    L.metrique = {};
    L.valeur = [];
    L.unite = {};
    % 1. Demarrage : [0 ; k_dem)
    k_dem = min([round(0.05 / Te), evts]);
    i10 = find(v >= 10, 1);
    i90 = find(v >= 90, 1);
    if isempty(i10) || isempty(i90)
        t_montee = NaN;
    else
        t_montee = (i90 - i10) * Te * 1e3;
    end
    L = ajouter(L, 'demarrage', 'demarrage', 't_montee', t_montee, 'ms');
    L = ajouter(L, 'demarrage', 'demarrage', 'depassement', max(0, (max(v(1:k_dem)) - 100) / 100 * 100), '%');
    hors = find(abs(e(1:k_dem)) > 1);
    if isempty(hors)
        t_etab = 0;
    else
        t_etab = hors(end) * Te * 1e3;                   % (k + 1) Te, k dernier instant hors bande
    end
    L = ajouter(L, 'demarrage', 'demarrage', 't_etablissement', t_etab, 'ms');
    L = ajouter(L, 'demarrage', 'demarrage', 'iL_crete', max(iL(1:k_dem)), 'A');
    % 2. Fenetre de classement [ka ; kb)
    if strcmp(code, 'S1')
        ka = 0;
        kb = k30;
    elseif strcmp(code, 'S10')
        ka = k100;
        kb = N;
    else
        ka = k30;
        kb = N;
    end
    idx = (ka + 1):kb;
    ew = e(idx);
    tau = ((ka:kb - 1)' - ka) * Te;                       % temps depuis le debut de la fenetre
    L = ajouter(L, 'classement', 'classement', 'IAE', sum(abs(ew)) * Te * 1e3, 'mV.s');
    L = ajouter(L, 'classement', 'classement', 'ISE', sum(ew .^ 2) * Te * 1e3, 'V2.ms');
    L = ajouter(L, 'classement', 'classement', 'ITAE', sum(tau .* abs(ew)) * Te * 1e3, 'mV.s2');
    L = ajouter(L, 'classement', 'classement', 'e_max', max(abs(ew)), 'V');
    L = ajouter(L, 'classement', 'classement', 'e_eff', sqrt(mean(ew .^ 2)) * 1e3, 'mV');
    % 4. Evenements (calcules ici pour compter ceux de la fenetre de
    % classement qui ne reviennent pas ; ajoutes a la liste plus bas)
    bornes = [evts, N];
    ev_iae = zeros(1, numel(evts));
    ev_emax = zeros(1, numel(evts));
    ev_ret = zeros(1, numel(evts));
    ev_rev = false(1, numel(evts));
    for j = 1:numel(evts)
        a = evts(j);
        b = bornes(j + 1);
        ef = e((a + 1):b);
        ev_iae(j) = sum(abs(ef)) * Te * 1e3;
        ev_emax(j) = max(abs(ef));
        hors = find(abs(ef) > 1);                         % positions p = k - a + 1
        ev_rev(j) = max(abs(e((max(a, b - n10) + 1):b))) <= 1;
        if isempty(hors)
            ev_ret(j) = 0;
        elseif ev_rev(j)
            ev_ret(j) = hors(end) * Te * 1e3;             % (k - a + 1) Te
        else
            ev_ret(j) = NaN;
        end
    end
    n_non = 0;
    for j = 1:numel(evts)
        if evts(j) >= ka && evts(j) < kb && ~ev_rev(j)
            n_non = n_non + 1;
        end
    end
    L = ajouter(L, 'classement', 'classement', 'n_non_revenus', n_non, '-');
    % 3. Regime permanent
    prec = [0, evts];
    for j = 1:numel(evts)
        a = max(prec(j), evts(j) - n5);
        b = evts(j);
        f = sprintf('avant_%gms', evts(j) * Te * 1e3);
        em = mean(e((a + 1):b));
        L = ajouter(L, 'regime_permanent', f, 'e_moy', em * 1e3, 'mV');
        L = ajouter(L, 'regime_permanent', f, 'e_moy_abs', abs(em) * 1e3, 'mV');
        L = ajouter(L, 'regime_permanent', f, 'ondulation_Vout', (max(v((a + 1):b)) - min(v((a + 1):b))) * 1e3, 'mV');
    end
    a = N - n10;
    em = mean(e((a + 1):N));
    L = ajouter(L, 'regime_permanent', 'fin', 'e_moy', em * 1e3, 'mV');
    L = ajouter(L, 'regime_permanent', 'fin', 'e_moy_abs', abs(em) * 1e3, 'mV');
    L = ajouter(L, 'regime_permanent', 'fin', 'ondulation_Vout', (max(v((a + 1):N)) - min(v((a + 1):N))) * 1e3, 'mV');
    % 4. Evenements (suite)
    for j = 1:numel(evts)
        f = sprintf('ev_%gms', evts(j) * Te * 1e3);
        L = ajouter(L, 'evenement', f, 'IAE', ev_iae(j), 'mV.s');
        L = ajouter(L, 'evenement', f, 'e_max', ev_emax(j), 'V');
        L = ajouter(L, 'evenement', f, 't_retour', ev_ret(j), 'ms');
        L = ajouter(L, 'evenement', f, 'revenu', double(ev_rev(j)), '-');
    end
    % 5. Commande et conduction, fenetre de classement
    dw = d(idx);
    L = ajouter(L, 'commande', 'classement', 'dd_moyen', mean(abs(diff(dw))), '-');
    L = ajouter(L, 'commande', 'classement', 'butee', mean(dw <= 0.01 + 1e-12 | dw >= 0.99 - 1e-12) * 100, '%');
    L = ajouter(L, 'commande', 'classement', 'dcm', mean(iL(idx) < 0.01) * 100, '%');
end

function L = ajouter(L, groupe, fenetre, metrique, val, unite)
    L.groupe{end + 1} = groupe;
    L.fenetre{end + 1} = fenetre;
    L.metrique{end + 1} = metrique;
    L.valeur(end + 1) = val;
    L.unite{end + 1} = unite;
end

function E = lire_resultats(racine, dossiers, modele, suffixes)
    % Toutes les entrees des fichiers resultats_<modele><suffixe>.mat des
    % dossiers donnes (cellule de structures) ; suffixes {'', '_S10'} pour
    % Simuler_<M>.m, {'_comparaison'} pour Simuler_Modeles_Comparaison.m.
    E = {};
    for i = 1:numel(dossiers)
        if ~isfolder(fullfile(racine, dossiers{i}))
            fprintf('  (dossier introuvable : %s)\n', fullfile(racine, dossiers{i}));
            continue;
        end
        for s = 1:numel(suffixes)
            f = fullfile(racine, dossiers{i}, ['resultats_' modele suffixes{s} '.mat']);
            if ~isfile(f)
                fprintf('  (absent : %s)\n', f);
                continue;
            end
            try
                S = load(f, 'resultats');
            catch ME
                fprintf('  ATTENTION : lecture impossible de %s : %s\n', f, ME.message);
                continue;
            end
            if ~isfield(S, 'resultats')
                fprintf('  ATTENTION : pas de variable "resultats" dans %s\n', f);
                continue;
            end
            for j = 1:numel(S.resultats)
                en.code = S.resultats(j).code;
                en.r = S.resultats(j);
                en.fichier = f;
                en.dossier = fullfile(racine, dossiers{i});
                en.priorite = strcmp(suffixes{s}, '_S10');   % fichier _S10 : prioritaire pour S10
                E{end + 1} = en; %#ok<AGROW>
            end
        end
    end
end

function [E, refus] = controler_correspondance(E, methode, BANC, seuil, tol_S10, sans_controle, Te)
    % Refuse les fichiers dont un scenario s'ecarte du banc (IAE de
    % classement, metriques_banc.csv) de plus du seuil : ce n'est pas la
    % methode retenue. S10 : seuil = max(seuil, tolerance S10
    % documentee). Scenarios de sans_controle et scenarios absents du banc :
    % non controles. Rend E sans les entrees des fichiers refuses et la liste
    % des refus (texte).
    refus = {};
    mauvais = {};
    for i = 1:numel(E)
        en = E{i};
        if any(strcmp(en.code, sans_controle)) || any(strcmp(mauvais, en.fichier))
            continue;
        end
        k = cle(en.code, methode, 'classement', 'classement', 'IAE');
        if ~isKey(BANC, k)
            continue;
        end
        s_ = seuil;
        if strcmp(en.code, 'S10') && ~isnan(tol_S10)
            s_ = max(seuil, tol_S10);
        end
        b = BANC(k);
        r = en.r;
        iae = NaN;
        if isequal(numel(r.v), numel(r.d), numel(r.iL))
            L = calculer_metriques(r.v, r.consigne, r.d, r.iL, lire_evenements(en), en.code, Te);
            q = find(strcmp(L.groupe, 'classement') & strcmp(L.metrique, 'IAE'), 1);
            iae = L.valeur(q);
        end
        if ~(abs(iae / b - 1) <= s_)
            mauvais{end + 1} = en.fichier; %#ok<AGROW>
            txt = sprintf(['%s : %s, ce fichier ne correspond pas a la methode du banc (IAE de classement %s ' ...
                           '%.2f mV.s contre %.2f au banc, %+.1f %%, seuil %.0f %%)'], methode, en.fichier, ...
                          en.code, iae, b, (iae / b - 1) * 100, s_ * 100);
            refus{end + 1} = txt; %#ok<AGROW>
            fprintf('  REFUSE : %s\n', txt);
        end
    end
    if ~isempty(mauvais)
        garde = true(1, numel(E));
        for i = 1:numel(E)
            garde(i) = ~any(strcmp(mauvais, E{i}.fichier));
        end
        E = E(garde);
        fprintf(['  -> %s : verifier le dossier lu (reglages DOSSIER_* en tete du script) ; ' ...
                 'pour l''ELM-PID, DOSSIER_ELM doit etre le dossier ELM_PID.\n'], methode);
    end
end

function en = choisir(E, code)
    % Entree de l'essai demande : la premiere trouvee, sauf si un fichier
    % _S10 en contient une aussi (prioritaire).
    en = [];
    for i = 1:numel(E)
        if strcmp(E{i}.code, code)
            if isempty(en) || (E{i}.priorite && ~en.priorite)
                en = E{i};
            end
        end
    end
end

function evts = lire_evenements(en)
    % Instants des evenements (s) : scenario_<code>.mat du dossier, sinon
    % les instants enregistres dans les grandeurs de Simuler_<M>.m.
    f = fullfile(en.dossier, ['scenario_' en.code '.mat']);
    if isfile(f)
        m = load(f, 'evenements');
        evts = reshape(double(m.evenements), 1, []);
    else
        g = en.r.grandeurs;
        if isempty(g.evenements)
            evts = zeros(1, 0);
        else
            evts = [g.evenements.t_ms] / 1e3;
        end
        fprintf('  (%s absent : evenements lus dans les grandeurs)\n', f);
    end
end

function k = cle(s, m, g, f, me)
    k = [s '|' m '|' g '|' f '|' me];
end

function x = valeur(BASE, s, m, g, f, me)
    k = cle(s, m, g, f, me);
    if isKey(BASE, k)
        x = BASE(k);
    else
        x = NaN;
    end
end

function t = nombre_csv(x)
    if isnan(x)
        t = 'NaN';
    else
        t = sprintf('%.10g', x);
    end
end

function B = lire_csv(f)
    % metriques_banc.csv -> table cle -> valeur.
    fid = fopen(f, 'r', 'n', 'UTF-8');
    C = textscan(fid, '%s %s %s %s %s %s %s', 'Delimiter', ';', 'Whitespace', '', 'HeaderLines', 1);
    fclose(fid);
    B = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:numel(C{1})
        B(cle(C{1}{i}, C{2}{i}, C{3}{i}, C{4}{i}, C{5}{i})) = str2double(C{6}{i});
    end
end

function T = ligne_md(T, cellules)
    T{end + 1} = ['| ' strjoin(cellules, ' | ') ' |'];
end

function t = f2(x, fmt)
    if isnan(x)
        t = 'n.d.';
    else
        t = sprintf(fmt, x);
    end
end

function txt = texte_classement(noms, vals, etoiles)
    % Ordre croissant, '=' entre voisines a moins de 5 %, '<' sinon ; '*' :
    % un evenement de la fenetre de classement n'est pas revenu.
    garde = ~isnan(vals);
    noms = noms(garde);
    vals = vals(garde);
    etoiles = etoiles(garde);
    if isempty(vals)
        txt = 'n.d.';
        return;
    end
    [vs, o] = sort(vals);
    txt = '';
    for i = 1:numel(o)
        if i > 1
            if vs(i) / vs(i - 1) - 1 < 0.05
                txt = [txt ' = ']; %#ok<AGROW>
            else
                txt = [txt ' < ']; %#ok<AGROW>
            end
        end
        txt = [txt sprintf('%s %.2f%s', noms{o(i)}, vs(i), etoiles{o(i)})]; %#ok<AGROW>
    end
end

function T = tableaux_communs(T, BASE, SCENARIOS, METHODES)
    % Tableaux identiques a ceux de metriques_banc.py (fonction
    % tableaux_communs).
    FEN = struct('S1', 'IAE de 0 a 30 ms', 'S2', 'IAE de 30 ms a la fin', 'S3', 'IAE de 30 ms a la fin', ...
                 'S8a', 'IAE de 30 ms a la fin', 'S10', 'IAE de 100 ms a la fin');
    % Resume
    T{end + 1} = '';
    T{end + 1} = '## Resume : IAE de classement (mV.s) ; erreur moyenne absolue en fin d''essai (mV)';
    T{end + 1} = '';
    T = ligne_md(T, [{'Methode'}, SCENARIOS]);
    T = ligne_md(T, repmat({'---'}, 1, numel(SCENARIOS) + 1));
    for im = 1:numel(METHODES)
        c = {METHODES{im}};
        for is = 1:numel(SCENARIOS)
            iae = valeur(BASE, SCENARIOS{is}, METHODES{im}, 'classement', 'classement', 'IAE');
            erp = valeur(BASE, SCENARIOS{is}, METHODES{im}, 'regime_permanent', 'fin', 'e_moy_abs');
            c{end + 1} = [f2(iae, '%.2f') ' ; ' f2(erp, '%.2f')]; %#ok<AGROW>
        end
        T = ligne_md(T, c);
    end
    % Classements
    T{end + 1} = '';
    T{end + 1} = '## Classement de chaque scenario (IAE de la fenetre de classement, mV.s ; "=" : moins de 5 % ; "*" : un evenement de la fenetre non revenu)';
    T{end + 1} = '';
    for is = 1:numel(SCENARIOS)
        s = SCENARIOS{is};
        if strcmp(s, 'S10')
            classees = METHODES(2:end);
        else
            classees = METHODES;
        end
        vals = zeros(1, numel(classees));
        etoiles = cell(1, numel(classees));
        for im = 1:numel(classees)
            vals(im) = valeur(BASE, s, classees{im}, 'classement', 'classement', 'IAE');
            nn = valeur(BASE, s, classees{im}, 'classement', 'classement', 'n_non_revenus');
            if ~isnan(nn) && nn > 0
                etoiles{im} = '*';
            else
                etoiles{im} = '';
            end
        end
        T{end + 1} = sprintf('- %s (%s) : %s', s, FEN.(s), texte_classement(classees, vals, etoiles)); %#ok<AGROW>
        if strcmp(s, 'S10')
            zn = valeur(BASE, s, METHODES{1}, 'classement', 'classement', 'IAE');
            nn = valeur(BASE, s, METHODES{1}, 'classement', 'classement', 'n_non_revenus');
            r = {};
            for im = 2:numel(METHODES)
                r{end + 1} = sprintf('%s %s', METHODES{im}, f2(valeur(BASE, s, METHODES{im}, 'classement', 'classement', 'IAE') / zn, '%.3f')); %#ok<AGROW>
            end
            T{end + 1} = sprintf('  Ziegler-Nichols (reference, hors classement) %s mV.s, %s evenement(s) de la fenetre non revenu(s) ; IAE / IAE de ZN : %s', ...
                                 f2(zn, '%.2f'), f2(nn, '%.0f'), strjoin(r, ', ')); %#ok<AGROW>
        end
    end
    % Par scenario : demarrage, fenetre de classement, commande
    for is = 1:numel(SCENARIOS)
        s = SCENARIOS{is};
        T{end + 1} = ''; %#ok<AGROW>
        T{end + 1} = sprintf('## %s : demarrage, fenetre de classement (%s), commande', s, FEN.(s)); %#ok<AGROW>
        T{end + 1} = ''; %#ok<AGROW>
        T = ligne_md(T, {'Methode', 'montee 10-90 % (ms)', 'depassement (%)', 'etabli +-1 V (ms)', 'iL crete (A)', ...
                         'IAE (mV.s)', 'ISE (V2.ms)', 'ITAE (mV.s2)', 'e max (V)', 'e eff (mV)', '|dd| moyen', ...
                         'butee (%)', 'DCM (%)'});
        T = ligne_md(T, repmat({'---'}, 1, 13));
        for im = 1:numel(METHODES)
            m = METHODES{im};
            T = ligne_md(T, {m, ...
                f2(valeur(BASE, s, m, 'demarrage', 'demarrage', 't_montee'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'demarrage', 'demarrage', 'depassement'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'demarrage', 'demarrage', 't_etablissement'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'demarrage', 'demarrage', 'iL_crete'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'classement', 'classement', 'IAE'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'classement', 'classement', 'ISE'), '%.4g'), ...
                f2(valeur(BASE, s, m, 'classement', 'classement', 'ITAE'), '%.4g'), ...
                f2(valeur(BASE, s, m, 'classement', 'classement', 'e_max'), '%.2f'), ...
                f2(valeur(BASE, s, m, 'classement', 'classement', 'e_eff'), '%.1f'), ...
                f2(valeur(BASE, s, m, 'commande', 'classement', 'dd_moyen'), '%.4f'), ...
                f2(valeur(BASE, s, m, 'commande', 'classement', 'butee'), '%.1f'), ...
                f2(valeur(BASE, s, m, 'commande', 'classement', 'dcm'), '%.1f')});
        end
    end
    % Regime permanent et evenements : fenetres lues dans la base
    for is = 1:numel(SCENARIOS)
        s = SCENARIOS{is};
        [fen_rp, fen_ev] = fenetres(BASE, s, METHODES);
        T{end + 1} = ''; %#ok<AGROW>
        T{end + 1} = sprintf('## %s : erreur en regime permanent, moyenne signee / moyenne absolue / ondulation crete a crete de Vout (mV)', s); %#ok<AGROW>
        T{end + 1} = ''; %#ok<AGROW>
        T = ligne_md(T, [{'Fenetre'}, METHODES]);
        T = ligne_md(T, repmat({'---'}, 1, numel(METHODES) + 1));
        for i = 1:numel(fen_rp)
            c = {texte_fenetre_rp(fen_rp{i})};
            for im = 1:numel(METHODES)
                c{end + 1} = sprintf('%s / %s / %s', ...
                    f2(valeur(BASE, s, METHODES{im}, 'regime_permanent', fen_rp{i}, 'e_moy'), '%+.2f'), ...
                    f2(valeur(BASE, s, METHODES{im}, 'regime_permanent', fen_rp{i}, 'e_moy_abs'), '%.2f'), ...
                    f2(valeur(BASE, s, METHODES{im}, 'regime_permanent', fen_rp{i}, 'ondulation_Vout'), '%.1f')); %#ok<AGROW>
            end
            T = ligne_md(T, c);
        end
        if isempty(fen_ev)
            continue;
        end
        T{end + 1} = ''; %#ok<AGROW>
        T{end + 1} = sprintf('## %s : par evenement, IAE (mV.s) / ecart max (V) / retour dans +-1 V', s); %#ok<AGROW>
        T{end + 1} = ''; %#ok<AGROW>
        T = ligne_md(T, [{'Evenement'}, METHODES]);
        T = ligne_md(T, repmat({'---'}, 1, numel(METHODES) + 1));
        for i = 1:numel(fen_ev)
            c = {strrep(strrep(fen_ev{i}, 'ev_', ''), 'ms', ' ms')};
            for im = 1:numel(METHODES)
                m = METHODES{im};
                iae = valeur(BASE, s, m, 'evenement', fen_ev{i}, 'IAE');
                if isnan(iae)
                    c{end + 1} = 'n.d.'; %#ok<AGROW>
                    continue;
                end
                ret = valeur(BASE, s, m, 'evenement', fen_ev{i}, 't_retour');
                if valeur(BASE, s, m, 'evenement', fen_ev{i}, 'revenu') == 0
                    tr = 'pas revenu';
                elseif ret == 0
                    tr = 'dans la bande';
                else
                    tr = sprintf('%.2f ms', ret);
                end
                c{end + 1} = sprintf('%.2f / %.2f / %s', iae, ...
                    valeur(BASE, s, m, 'evenement', fen_ev{i}, 'e_max'), tr); %#ok<AGROW>
            end
            T = ligne_md(T, c);
        end
    end
end

function [fen_rp, fen_ev] = fenetres(BASE, s, METHODES)
    % Noms des fenetres de regime permanent et des evenements d'un scenario,
    % dans l'ordre des instants, pris chez la premiere methode presente.
    fen_rp = {};
    fen_ev = {};
    cles = keys(BASE);
    for im = 1:numel(METHODES)
        pre_rp = [s '|' METHODES{im} '|regime_permanent|'];
        pre_ev = [s '|' METHODES{im} '|evenement|'];
        for i = 1:numel(cles)
            k = cles{i};
            if startsWith(k, pre_rp) && endsWith(k, '|e_moy')
                fen_rp{end + 1} = k(numel(pre_rp) + 1:end - numel('|e_moy')); %#ok<AGROW>
            elseif startsWith(k, pre_ev) && endsWith(k, '|IAE')
                fen_ev{end + 1} = k(numel(pre_ev) + 1:end - numel('|IAE')); %#ok<AGROW>
            end
        end
        if ~isempty(fen_rp)
            break;
        end
    end
    fen_rp = trier_fenetres(fen_rp);
    fen_ev = trier_fenetres(fen_ev);
end

function f = trier_fenetres(f)
    % Tri par instant ('avant_115ms', 'ev_115ms') ; 'fin' a la fin.
    t = zeros(1, numel(f));
    for i = 1:numel(f)
        if strcmp(f{i}, 'fin')
            t(i) = inf;
        else
            t(i) = str2double(regexprep(f{i}, '^(avant_|ev_)|ms$', ''));
        end
    end
    [~, o] = sort(t);
    f = f(o);
end

function t = texte_fenetre_rp(f)
    if strcmp(f, 'fin')
        t = '10 dernieres ms';
    else
        t = ['5 ms ' strrep(strrep(f, 'avant_', 'avant '), 'ms', ' ms')];
    end
end
