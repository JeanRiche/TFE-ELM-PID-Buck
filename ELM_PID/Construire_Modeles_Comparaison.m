%% CONSTRUIRE_MODELES_COMPARAISON.m
%
% OBJECTIF
% --------
% Construire, pour les scenarios du corps du memoire (S1, S2, S3, S8a,
% S10), un modele Simulink par scenario qu'on ouvre et qu'on lance d'un
% clic (bouton Run), sans script prealable :
%     <PREFIXE>_<code>.slx     par exemple ELM_PID_S10.slx, ZN_S8a.slx
% Le script est le meme dans les quatre dossiers de methode (ELM_PID,
% PINN_PID, PSO_PID, FUZZY_PID). Il trouve seul le modele de la methode du
% dossier : le seul fichier Buck_Commun_*.slx (Buck_Commun_ELM_PID.slx,
% Buck_Commun_PINN_PID.slx, Buck_Commun_PSO_PID.slx ou
% Buck_Commun_Fuzzy_PID.slx). Avec AVEC_ZIEGLER_NICHOLS = true, il
% construit aussi les modeles de Ziegler-Nichols (ZN_<code>.slx) a partir
% de Buck_Commun.slx ; un seul dossier suffit pour eux.
%
% CE QUE LE SCRIPT CHANGE, ET RIEN D'AUTRE
% -----------------------------------------
% Chaque modele est une copie fraiche du modele source (fichier copie,
% puis charge sous son nouveau nom). Seuls changent :
%   1. le rappel InitFcn : il se place dans le dossier du modele, puis
%      charge le scenario sans condition, charger_scenario('<code>') ;
%      (le rappel du modele commun, qui charge S1 si aucun scenario n'est
%      en memoire, est remplace ; un rappel different, inattendu, est garde
%      a la suite) ;
%   2. le StopTime : un demi-pas Tc apres le dernier instant du scenario,
%      (n - 0.5) x Tc avec n = nombre d'instants du profil, comme le
%      StopTime que Simuler_<Methode>.m passe a sim() ;
%   3. des blocs d'observation, qui ne renvoient rien dans le modele :
%        Scope "Observation Vout et consigne"   (Vout et 100 V + dvref)
%        Scope "Observation iL"                 (courant de la bobine)
%        Scope "Observation rapport cyclique"   (signal qui entre dans le PWM)
%        Gain  "Observation gains rapportes a ZN" puis
%        Scope "Observation gains"              (methodes adaptatives :
%                                                P, I, D divises par ceux
%                                                de Ziegler-Nichols)
%        To Workspace "Enregistrement Vout", "Enregistrement iL",
%        "Enregistrement rapport cyclique" et (methodes adaptatives)
%        "Enregistrement gains" : variables cmp_vout, cmp_iL, cmp_d,
%        cmp_K, echantillonnees a Tc = 1/(22000*10) s (texte exact),
%        format Timeseries, sans limite de points.
%      Les Scope sont echantillonnes a Tc eux aussi (44 001 points par
%      essai de 0.2 s au lieu de 5.3 millions au pas du powergui).
% Le circuit, le regulateur, ses gains, le PWM, le powergui et le reste
% de la configuration ne changent pas. Les deux To Workspace d'origine du
% modele commun (variables x et y, au pas du powergui) sont gardes tels
% quels.
%
% APRES LA CONSTRUCTION
% ---------------------
% Ouvrir par exemple ELM_PID_S10.slx (double-clic dans le dossier courant
% de MATLAB) et cliquer sur Run. Le rappel InitFcn charge S10 ; a la fin,
% cmp_vout, cmp_iL, cmp_d (et cmp_K) sont dans le workspace, et les Scope
% montrent les courbes. Le temps de calcul et les resultats au format de
% Simuler_<Methode>.m : Simuler_Modeles_Comparaison.m ; les metriques
% (definitions de COMPARAISON/metriques/metriques_banc.py) :
% COMPARAISON/metriques/Metriques_Simulink.m avec RESULTATS = 'comparaison'.
%
% Prerequis dans le dossier courant MATLAB : le modele de la methode
% (Buck_Commun_<Methode>.slx, construit par Construction_<Methode>.m),
% Buck_Commun.slx (pour Ziegler-Nichols), charger_scenario.m,
% scenario_S1.mat, scenario_S2.mat, scenario_S3.mat, scenario_S8a.mat,
% scenario_S10.mat et les fichiers du bloc de la methode (par exemple
% elm_pid_adaptatif.m, elm_pid_modele.mat, elm_pid_reglages.mat,
% ensemble_gains_elm.mat pour l'ELM-PID).
% Duree : moins d'une minute. Compatible MATLAB R2024a.
% Les scripts valides (Construction_*.m, Simuler_*.m) ne sont pas touches,
% et les modeles sources ne sont jamais modifies.
% ---------------------------------------------------------------------

clear; clc;

% --- Reglages ---
CODES = {'S1', 'S2', 'S3', 'S8a', 'S10'};        % scenarios du corps du memoire
AVEC_METHODE = true;                             % modeles de la methode du dossier
AVEC_ZIEGLER_NICHOLS = false;                    % true : aussi ZN_<code>.slx (depuis Buck_Commun.slx)

TC_TXT = '1/(22000*10)';                         % periode du regulateur (texte exact)
TC = 1/(22000*10);                               % la meme, en nombre (s)
K_ZN = [0.093910, 301.089, 7.3227e-06];          % P, I, D de Ziegler-Nichols (pour la vue des gains)
PWM = sprintf('PWM Generator\n(DC-DC)');         % nom du bloc PWM (vrai saut de ligne)
INIT_BASE = 'if ~exist(''sc_R0'', ''var''), charger_scenario(''S1''); end';   % rappel de Buck_Commun.slx

% Methodes connues : modele source, prefixe des modeles construits, bloc dont
% la sortie 1 est le vecteur des gains [P I D] ('' : gains constants), et
% fichiers dont le bloc a besoin.
METHODES = { ...
    'Buck_Commun_ELM_PID',   'ELM_PID',   'ELM-PID Adaptatif', ...
        {'elm_pid_adaptatif.m', 'elm_pid_modele.mat', 'elm_pid_reglages.mat', 'ensemble_gains_elm.mat'}; ...
    'Buck_Commun_PINN_PID',  'PINN_PID',  'PINN-PID Adaptatif', ...
        {'pinn_pid_adaptatif.m', 'pinn_pid_modele.mat', 'pinn_pid_reglages.mat'}; ...
    'Buck_Commun_Fuzzy_PID', 'Fuzzy_PID', 'Ordonnanceur Flou Zhao', ...
        {'ordonnanceur_flou_zhao.m', 'fuzzy_pid_reglages.mat'}; ...
    'Buck_Commun_PSO_PID',   'PSO_PID',   '', {}; ...
    'Buck_Commun',           'ZN',        '', {}};

% --- Modeles a construire ---
travaux = {};                                    % lignes de METHODES retenues
if AVEC_METHODE
    trouves = dir(fullfile(pwd, 'Buck_Commun_*.slx'));
    noms = erase({trouves.name}, '.slx');
    if isempty(noms) && ~AVEC_ZIEGLER_NICHOLS
        error(['Aucun modele de methode (Buck_Commun_*.slx) dans le dossier courant %s.\n' ...
               'Lancer d''abord Construction_<Methode>.m, ou mettre AVEC_ZIEGLER_NICHOLS = true ' ...
               'pour ne construire que les modeles de Ziegler-Nichols.'], pwd);
    elseif numel(noms) > 1
        error(['Plusieurs modeles de methode dans le dossier courant : %s.\n' ...
               'Il doit y en avoir un seul (le script est fait pour un dossier de methode).'], strjoin(noms, ', '));
    elseif numel(noms) == 1
        i = find(strcmp(METHODES(:, 1), noms{1}));
        if isempty(i)
            error(['%s.slx : methode inconnue de ce script. Methodes connues : %s.\n' ...
                   'Ajouter une ligne au tableau METHODES en tete du script.'], noms{1}, ...
                  strjoin(METHODES(1:end-1, 1)', ', '));
        end
        travaux{end + 1} = i;
        fprintf('Modele de la methode du dossier : %s.slx (modeles %s_<code>.slx).\n', noms{1}, METHODES{i, 2});
    else
        fprintf('Pas de modele de methode dans ce dossier : seuls les modeles de Ziegler-Nichols seront construits.\n');
    end
end
if AVEC_ZIEGLER_NICHOLS
    travaux{end + 1} = size(METHODES, 1);
end
if isempty(travaux)
    error('Rien a construire : AVEC_METHODE et AVEC_ZIEGLER_NICHOLS sont tous deux a false.');
end

% --- Prerequis communs ---
manquants = {};
for f = [{'charger_scenario.m'}, strcat('scenario_', CODES, '.mat')]
    if ~isfile(fullfile(pwd, f{1}))
        manquants{end + 1} = f{1}; %#ok<SAGROW>
    end
end
for w = travaux
    ligne = METHODES(w{1}, :);
    for f = [{[ligne{1} '.slx']}, ligne{4}]
        if ~isfile(fullfile(pwd, f{1}))
            manquants{end + 1} = f{1}; %#ok<SAGROW>
        end
    end
end
if ~isempty(manquants)
    error(['Fichier(s) introuvable(s) dans le dossier courant %s :\n  %s\n' ...
           'Se placer dans le dossier de la methode (Current Folder) ; les scenario_*.mat et ' ...
           'charger_scenario.m sont fournis avec la base commune, les fichiers du bloc avec la methode.'], ...
          pwd, strjoin(manquants, '\n  '));
end
for f = {'charger_scenario.m'}
    copies = which(f{1}, '-all');
    if numel(copies) ~= 1
        error('Il doit exister exactement un %s sur le chemin MATLAB, il y en a %d :\n%s', ...
              f{1}, numel(copies), strjoin(copies, newline));
    end
end

% --- Construction ---
construits = {};
for w = travaux
    ligne = METHODES(w{1}, :);
    SOURCE = ligne{1};
    PREFIXE = ligne{2};
    BLOC_GAINS = ligne{3};
    if bdIsLoaded(SOURCE) && strcmp(get_param(SOURCE, 'Dirty'), 'on')
        error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', SOURCE);
    end
    copies = which([SOURCE '.slx'], '-all');
    if numel(copies) ~= 1
        error('Il doit exister exactement un %s.slx sur le chemin MATLAB, il y en a %d :\n%s', ...
              SOURCE, numel(copies), strjoin(copies, newline));
    end
    fprintf('\n=== %s.slx -> %s_<code>.slx ===\n', SOURCE, PREFIXE);
    for c = CODES
        code = c{1};
        MDL = [PREFIXE '_' code];
        chemin = fullfile(pwd, [MDL '.slx']);
        sc = charger_scenario(code);                     % verifie aussi la periode des profils
        n = numel(sc.t);
        stop_txt = sprintf('%.17g', (n - 0.5) * TC);     % meme StopTime que Simuler_<Methode>.m
        if bdIsLoaded(MDL)
            if strcmp(get_param(MDL, 'Dirty'), 'on')
                error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', MDL);
            end
            close_system(MDL, 0);
        end
        if isfile(chemin)
            delete(chemin);
            fprintf('  ancien %s.slx supprime (on repart toujours d''une copie fraiche)\n', MDL);
        end
        copyfile(fullfile(pwd, [SOURCE '.slx']), chemin);
        fileattrib(chemin, '+w');                        % la copie doit pouvoir etre enregistree
        copies = which([MDL '.slx'], '-all');
        if numel(copies) ~= 1
            delete(chemin);
            error('Il doit exister exactement un %s.slx sur le chemin MATLAB, il y en a %d :\n%s', ...
                  MDL, numel(copies), strjoin(copies, newline));
        end
        try
            load_system(chemin);
            % Blocs observes (ils doivent exister dans le modele source)
            vout = [MDL '/Vout'];
            il = [MDL '/iL'];
            consigne = [MDL '/Somme Consigne'];
            for b = {vout, il, consigne, [MDL '/' PWM]}
                if getSimulinkBlockHandle(b{1}) == -1
                    error('Bloc introuvable dans %s.slx : %s', SOURCE, strrep(b{1}, newline, ' '));
                end
            end
            if ~isempty(BLOC_GAINS) && getSimulinkBlockHandle([MDL '/' BLOC_GAINS]) == -1
                error('Bloc des gains introuvable dans %s.slx : %s', SOURCE, BLOC_GAINS);
            end
            src_pwm = source_entree_pwm(MDL, PWM);       % ce qui alimente l'entree du PWM

            % 1. Rappel InitFcn
            ancien = get_param(MDL, 'InitFcn');
            set_param(MDL, 'InitFcn', texte_initfcn(MDL, code, ancien, INIT_BASE));
            if ~strcmp(strtrim(ancien), INIT_BASE) && ~isempty(strtrim(ancien))
                warning('%s : rappel InitFcn inattendu dans %s.slx, garde a la suite du chargement du scenario :\n%s', ...
                        MDL, SOURCE, ancien);
            end
            % 2. StopTime
            set_param(MDL, 'StopTime', stop_txt);
            % 3. Blocs d'observation, a droite du schema
            blocs = find_system(MDL, 'SearchDepth', 1, 'Type', 'Block');
            pos = cell2mat(get_param(blocs, 'Position'));
            x0 = max(pos(:, 3)) + 250;
            y0 = min(pos(:, 2));
            ajouter_scope(MDL, 'Observation Vout et consigne', [x0, y0, x0 + 60, y0 + 70], ...
                          {vout, 1; consigne, 1}, TC_TXT);
            ajouter_to_workspace(MDL, 'Enregistrement Vout', 'cmp_vout', [x0 + 140, y0 + 20, x0 + 230, y0 + 50], ...
                                 vout, 1, TC_TXT);
            ajouter_scope(MDL, 'Observation iL', [x0, y0 + 110, x0 + 60, y0 + 170], {il, 1}, TC_TXT);
            ajouter_to_workspace(MDL, 'Enregistrement iL', 'cmp_iL', [x0 + 140, y0 + 125, x0 + 230, y0 + 155], ...
                                 il, 1, TC_TXT);
            ajouter_scope(MDL, 'Observation rapport cyclique', [x0, y0 + 210, x0 + 60, y0 + 270], ...
                          {src_pwm.bloc, src_pwm.port}, TC_TXT);
            ajouter_to_workspace(MDL, 'Enregistrement rapport cyclique', 'cmp_d', ...
                                 [x0 + 140, y0 + 225, x0 + 230, y0 + 255], src_pwm.bloc, src_pwm.port, TC_TXT);
            if ~isempty(BLOC_GAINS)
                gains = [MDL '/' BLOC_GAINS];
                gain_zn = [MDL '/Observation gains rapportes a ZN'];
                add_block('simulink/Math Operations/Gain', gain_zn, 'Position', [x0 - 110, y0 + 320, x0 - 50, y0 + 350]);
                set_param(gain_zn, 'Gain', sprintf('[1/%.17g, 1/%.17g, 1/%.17g]', K_ZN), ...
                          'Multiplication', 'Element-wise(K.*u)');
                relier(MDL, gains, 1, gain_zn, 1);
                ajouter_scope(MDL, 'Observation gains', [x0, y0 + 310, x0 + 60, y0 + 370], {gain_zn, 1}, TC_TXT);
                ajouter_to_workspace(MDL, 'Enregistrement gains', 'cmp_K', [x0 + 140, y0 + 325, x0 + 230, y0 + 355], ...
                                     gains, 1, TC_TXT);
            end
            save_system(MDL);
            close_system(MDL, 0);
        catch ME
            if bdIsLoaded(MDL)
                close_system(MDL, 0);                    % sans enregistrer
            end
            if isfile(chemin)
                delete(chemin);                          % pas de modele a moitie construit
            end
            rethrow(ME);
        end
        construits{end + 1} = MDL; %#ok<SAGROW>
        fprintf('  %s.slx : scenario %s (%s), StopTime %s s, %d blocs d''observation\n', MDL, code, sc.nom, ...
                stop_txt, 6 + 3 * ~isempty(BLOC_GAINS));
    end
end

fprintf(['\n%d modele(s) construit(s) : %s.\n' ...
         'Pour en lancer un : l''ouvrir (double-clic) et cliquer sur Run.\n' ...
         'Simulation et temps : Simuler_Modeles_Comparaison.m.\n'], ...
        numel(construits), strjoin(construits, ', '));


%% ===================== Fonctions =====================================

function txt = texte_initfcn(mdl, code, ancien, init_base)
    % Rappel InitFcn du modele de comparaison : se placer dans le dossier du
    % modele (charger_scenario lit scenario_<code>.mat dans le dossier
    % courant, et les blocs MATLAB System y lisent leurs fichiers), puis
    % charger le scenario, sans condition.
    lignes = {sprintf('%% Modele de comparaison %s (Construire_Modeles_Comparaison.m) : scenario %s.', mdl, code), ...
              sprintf('cmp_dossier_modele = fileparts(get_param(''%s'', ''FileName''));', mdl), ...
              'if ~strcmp(pwd, cmp_dossier_modele), cd(cmp_dossier_modele); end', ...
              'clear cmp_dossier_modele', ...
              sprintf('charger_scenario(''%s'');', code)};
    txt = strjoin(lignes, newline);
    if ~strcmp(strtrim(ancien), init_base) && ~isempty(strtrim(ancien))
        txt = [txt newline ancien];                      % rappel inattendu : garde a la suite
    end
end

function s = source_entree_pwm(mdl, pwm_nom)
    % Bloc et numero de port qui alimentent l'entree du PWM Generator.
    pwm = [mdl '/' pwm_nom];
    ph = get_param(pwm, 'PortHandles');
    ligne = get_param(ph.Inport(1), 'Line');
    if ligne == -1
        error('L''entree du PWM n''est reliee a rien dans %s.', mdl);
    end
    s.bloc = getfullname(get_param(ligne, 'SrcBlockHandle'));
    s.port = get_param(get_param(ligne, 'SrcPortHandle'), 'PortNumber');
end

function relier(mdl, bloc_src, port_src, bloc_dst, port_dst)
    % Ligne de signal entre deux ports (handles demandes au moment meme).
    ph_src = get_param(bloc_src, 'PortHandles');
    ph_dst = get_param(bloc_dst, 'PortHandles');
    add_line(mdl, ph_src.Outport(port_src), ph_dst.Inport(port_dst), 'autorouting', 'on');
end

function ajouter_scope(mdl, nom, position, sources, Te_txt)
    % Scope a autant d'entrees que de sources ({bloc, port; ...}), echantillonne
    % a Te_txt, toutes les courbes sur un meme axe.
    chemin = [mdl '/' nom];
    add_block('simulink/Sinks/Scope', chemin, 'Position', position);
    set_param(chemin, 'NumInputPorts', num2str(size(sources, 1)));
    try
        set_param(chemin, 'SampleTime', Te_txt);         % parametre du bloc Scope (R2024a)
    catch
        try
            cfg = get_param(chemin, 'ScopeConfiguration');
            cfg.SampleTime = Te_txt;                     % meme reglage par l'objet de configuration
        catch
            warning(['%s : periode d''echantillonnage du Scope non reglee (il heritera du pas du powergui ; ' ...
                     'la simulation reste la meme, l''affichage est seulement plus lourd).'], nom);
        end
    end
    try                                                  % presentation seulement
        cfg = get_param(chemin, 'ScopeConfiguration');
        cfg.LayoutDimensions = [1 1];
        cfg.ShowLegend = true;
        cfg.ShowGrid = true;
        cfg.Title = nom;
    catch
        % sans consequence sur la simulation
    end
    for i = 1:size(sources, 1)
        relier(mdl, sources{i, 1}, sources{i, 2}, chemin, i);
    end
end

function ajouter_to_workspace(mdl, nom, variable, position, bloc_src, port_src, Te_txt)
    % To Workspace (Timeseries, sans limite de points) echantillonne a Te_txt.
    chemin = [mdl '/' nom];
    add_block('simulink/Sinks/To Workspace', chemin, 'Position', position);
    set_param(chemin, 'VariableName', variable, 'SaveFormat', 'Timeseries', ...
              'SampleTime', Te_txt, 'MaxDataPoints', 'inf');
    relier(mdl, bloc_src, port_src, chemin, 1);
end
