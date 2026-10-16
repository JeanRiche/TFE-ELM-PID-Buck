%% SIMULER_MODELES_COMPARAISON.m
%
% OBJECTIF
% --------
% Simuler les modeles qu'on lance d'un clic, construits par
% Construire_Modeles_Comparaison.m (<PREFIXE>_<code>.slx, par exemple
% ELM_PID_S10.slx ou ZN_S8a.slx), mesurer le temps d'execution de chaque
% simulation, et ecrire les signaux au format des fichiers de
% Simuler_<Methode>.m pour que COMPARAISON/metriques/Metriques_Simulink.m
% calcule TOUTES les metriques avec les definitions de reference
% (en-tete de COMPARAISON/metriques/metriques_banc.py) :
% demarrage, IAE / ISE / ITAE de la fenetre de classement, ecarts, erreur
% en regime permanent (5 ms avant chaque evenement, 10 dernieres ms),
% ondulation, evenements, commande, temps de calcul Simulink.
% Ce script ne redefinit aucune metrique. Il affiche seulement, comme
% Simuler_<Methode>.m, les grandeurs du banc commun (fonction grandeurs,
% copie de Simuler_ELM_PID.m) et l'ecart au banc, pour reperer tout de
% suite un modele qui ne suit pas.
%
% Le meme script est dans les quatre dossiers de methode. Il trouve seul
% les modeles <PREFIXE>_<code>.slx du dossier courant.
%
% LE TEMPS D'EXECUTION
% --------------------
% duree_calcul_s = duree de sim() mesuree par tic / toc (compilation du
% modele et rappel InitFcn compris, comme un clic sur Run) ; us_par_pas =
% duree_calcul_s divisee par le nombre de periodes du regulateur (44 001
% pour 0.2 s). Memes noms de champs que Simuler_<Methode>.m : c'est le
% temps de calcul Simulink de Metriques_Simulink.m. Ce n'est pas le cout du
% regulateur : l'execution interpretee des blocs MATLAB System et le
% circuit au pas du powergui dominent.
% Pendant la simulation, les deux To Workspace d'origine du modele commun
% (variables x et y, au pas du powergui, 5.3 millions de points par bloc)
% sont mis en commentaire EN MEMOIRE, comme dans Simuler_<Methode>.m
% (reglage NEUTRALISER_TO_WORKSPACE_DU_MODELE) : les durees se comparent
% alors a celles de Simuler_<Methode>.m. Les blocs d'observation ajoutes
% (Scope, To Workspace cmp_*) restent actifs. Aucun .slx n'est enregistre.
%
% CE QUE PRODUIT CE SCRIPT
% -------------------------
%   - dans la console : par modele, le temps de calcul, les grandeurs du
%     banc commun et l'ecart au banc (si predictions_banc_*.json est la) ;
%     puis un tableau recapitulatif ;
%   - resultats_<MODELE SOURCE>_comparaison.mat, par exemple
%     resultats_Buck_Commun_ELM_PID_comparaison.mat ou, pour ZN,
%     resultats_Buck_Commun_comparaison.mat : variable "resultats", une
%     entree par scenario, memes champs que Simuler_<Methode>.m (code, nom,
%     t, v, consigne, mu, d, iL, K, grandeurs, duree_calcul_s, us_par_pas)
%     plus "modele" (nom du modele simule). Si le fichier existe deja, les
%     scenarios simules y sont remplaces et les autres gardes : on peut
%     simuler un scenario a la fois. Les fichiers de Simuler_<Methode>.m
%     (resultats_Buck_Commun_<M>.mat et _S10.mat) ne sont jamais touches.
%   - une figure par modele simule (reglage FIGURES) : quatre
%     tuiles empilees, abscisses en ms sur toute la duree, liees entre
%     elles ; lignes verticales pointillees aux evenements du scenario ;
%       1. Vout et consigne ; 2. rapport cyclique (signal qui entre dans le
%       PWM) ; 3. courant iL ; 4. gains P, I, D divises par ceux de
%       Ziegler-Nichols (methodes adaptatives, enregistrement cmp_K). Le
%       PSO-PID et Ziegler-Nichols n'ont pas de gains variables : la tuile 4
%       trace alors les gains constants lus dans le bloc "PID Controller"
%       divises par ceux de ZN, en droites horizontales pleines sans
%       marqueur (couleurs par defaut, P, I, D dans le meme ordre que les
%       methodes adaptatives, legende au-dessus). Chaque droite porte sa
%       valeur ecrite dans la tuile, par exemple "P = 0.23958 (2.551 x ZN)",
%       a 15, 45 ou 75 % de la duree, au-dessus ou au-dessous de la droite
%       de facon que les etiquettes de deux droites proches ne se
%       recouvrent pas. Rapports egaux a 5e-4 pres en relatif (ZN : tous
%       valent 1, meme avec des gains arrondis dans le bloc) : une seule
%       droite et une seule etiquette, "P, I, D = 1.000 x ZN (P = ...,
%       I = ..., D = ...)". Ordonnees elargies (marge de 35 % au moins)
%       pour que les etiquettes tiennent dans la tuile sans toucher les
%       bornes ; graduations aux valeurs des droites quand elles sont assez
%       espacees. Le texte des gains et de leur rapport a ceux de ZN reste
%       ecrit sous la tuile, avec l'axe des temps ; la figure est alors plus
%       haute (900 pixels au plus) et posee plus haut sur l'ecran pour que
%       ce texte ne soit pas cache Si les
%       gains ne se lisent pas, la tuile 4 reste un texte seul.
%     Titre : modele, methode, scenario. Echelle des ordonnees de chaque
%     tuile : du minimum au maximum des donnees tracees dans la tuile,
%     elargi de MARGE_ORDONNEES (8 %) de l'etendue de chaque cote ; aucune
%     courbe ne touche les bornes. Legendes au-dessus des tuiles, hors de la
%     zone des courbes.
%     Chaque figure est enregistree dans le sous-dossier figures_comparaison
%     du dossier courant : <modele>.png (exportgraphics, 200 dpi) et
%     <modele>.fig. FERMER_FIGURES = true ferme chaque figure apres
%     l'enregistrement (utile quand on simule beaucoup de modeles).
%     La figure est tracee APRES la mesure du temps de sim() : elle ne
%     change ni le temps de calcul, ni les grandeurs, ni la comparaison au
%     banc. Si le trace echoue, un avertissement est affiche et le script
%     continue (les resultats .mat sont ecrits quand meme).
% Ziegler-Nichols : les modeles ZN_<code>.slx (Construire_Modeles_Comparaison.m
% avec AVEC_ZIEGLER_NICHOLS = true, dans un seul dossier, par exemple
% PSO_PID) se simulent avec ce meme script, figures comprises :
%   AVEC_ZIEGLER_NICHOLS = true   -> modeles de la methode du dossier ET ZN ;
%   MODELES_A_SIMULER = {'ZN_S1', 'ZN_S2', 'ZN_S3', 'ZN_S8a', 'ZN_S10'}
%                                 -> Ziegler-Nichols seul.
% Resultats : resultats_Buck_Commun_comparaison.mat, figures ZN_<code>.png.
% Ensuite, quand les cinq methodes sont simulees : dossier courant
% COMPARAISON/metriques, Metriques_Simulink.m avec RESULTATS =
% 'comparaison' (metriques_simulink_comparaison.csv et .md).
%
% Prerequis dans le dossier courant MATLAB : les modeles
% <PREFIXE>_<code>.slx, charger_scenario.m, scenario_<code>.mat, les
% fichiers du bloc de la methode ; facultatif : predictions_banc_*.json.
% Duree par scenario de 0.2 s : environ 30 s pour ZN et le PSO-PID, une a
% deux minutes pour le Fuzzy-PID et l'ELM-PID, deux a quatre minutes pour
% le PINN-PID (execution interpretee des blocs MATLAB System).
% Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

clear; clc;

% --- Reglages ---
CODES = {'S1', 'S2', 'S3', 'S8a', 'S10'};        % scenarios cherches quand la liste est vide
MODELES_A_SIMULER = {};                          % {} : tous les <PREFIXE>_<code>.slx du dossier ;
                                                 % sinon par exemple {'ELM_PID_S10', 'ZN_S10'}
AVEC_ZIEGLER_NICHOLS = false;                    % liste vide : true pour simuler aussi les ZN_<code>.slx
NEUTRALISER_TO_WORKSPACE_DU_MODELE = true;       % To Workspace d'origine (x, y) en commentaire, en memoire
TC = 1/(22000*10);                               % periode du regulateur et des enregistrements (s)
FIGURES = true;                                  % une figure par modele simule (PNG et .fig enregistres)
FERMER_FIGURES = false;                          % true : fermer chaque figure une fois enregistree
DOSSIER_FIGURES = 'figures_comparaison';         % sous-dossier du dossier courant
MARGE_ORDONNEES = 0.08;                          % marge en ordonnee de chaque cote (part de l'etendue, >= 0.05)
K_ZN = [0.093910, 301.089, 7.3227e-06];          % P, I, D de Ziegler-Nichols (tuile des gains)

% Prefixes : prefixe des modeles, modele source (nom des fichiers de
% resultats, comme Simuler_<Methode>.m), methode, resultats attendus du banc.
PREFIXES = {'ELM_PID',   'Buck_Commun_ELM_PID',   'ELM-PID',         'predictions_banc_elm_pid.json'; ...
            'PINN_PID',  'Buck_Commun_PINN_PID',  'PINN-PID',        'predictions_banc_pinn_pid.json'; ...
            'Fuzzy_PID', 'Buck_Commun_Fuzzy_PID', 'Fuzzy-PID',       'predictions_banc_fuzzy_pid.json'; ...
            'PSO_PID',   'Buck_Commun_PSO_PID',   'PSO-PID',         'predictions_banc_pso_pid.json'; ...
            'ZN',        'Buck_Commun',           'Ziegler-Nichols', 'predictions_banc_pid_classique.json'};

% --- Liste des modeles ---
liste = struct('modele', {}, 'source', {}, 'methode', {}, 'code', {}, 'predictions', {});
if isempty(MODELES_A_SIMULER)
    for p = 1:size(PREFIXES, 1)
        if strcmp(PREFIXES{p, 1}, 'ZN') && ~AVEC_ZIEGLER_NICHOLS
            continue
        end
        for c = CODES
            mdl = [PREFIXES{p, 1} '_' c{1}];
            if isfile(fullfile(pwd, [mdl '.slx']))
                liste(end + 1) = struct('modele', mdl, 'source', PREFIXES{p, 2}, 'methode', PREFIXES{p, 3}, ...
                                        'code', c{1}, 'predictions', PREFIXES{p, 4}); %#ok<SAGROW>
            end
        end
    end
    if isempty(liste)
        error(['Aucun modele <PREFIXE>_<code>.slx dans le dossier courant %s.\n' ...
               'Lancer d''abord Construire_Modeles_Comparaison.m.'], pwd);
    end
else
    for m = reshape(cellstr(MODELES_A_SIMULER), 1, [])
        p = find(cellfun(@(x) startsWith(m{1}, [x '_']), PREFIXES(:, 1)), 1);
        if isempty(p)
            error('MODELES_A_SIMULER : %s ne commence par aucun prefixe connu (%s).', m{1}, ...
                  strjoin(PREFIXES(:, 1)', ', '));
        end
        if ~isfile(fullfile(pwd, [m{1} '.slx']))
            error('Modele introuvable dans le dossier courant : %s.slx (lancer Construire_Modeles_Comparaison.m).', m{1});
        end
        liste(end + 1) = struct('modele', m{1}, 'source', PREFIXES{p, 2}, 'methode', PREFIXES{p, 3}, ...
                                'code', m{1}(numel(PREFIXES{p, 1}) + 2:end), 'predictions', PREFIXES{p, 4}); %#ok<SAGROW>
    end
end
fprintf('Modeles simules : %s.\n', strjoin({liste.modele}, ', '));

% --- Prerequis ---
manquants = {};
for f = [{'charger_scenario.m'}, strcat('scenario_', unique({liste.code}, 'stable'), '.mat')]
    if ~isfile(fullfile(pwd, f{1}))
        manquants{end + 1} = f{1}; %#ok<SAGROW>
    end
end
if ~isempty(manquants)
    error('Fichier(s) introuvable(s) dans le dossier courant %s :\n  %s', pwd, strjoin(manquants, '\n  '));
end
for f = [{'charger_scenario.m'}, strcat({liste.modele}, '.slx')]
    copies = which(f{1}, '-all');
    if numel(copies) ~= 1
        error('Il doit exister exactement un %s sur le chemin MATLAB, il y en a %d :\n%s', ...
              f{1}, numel(copies), strjoin(copies, newline));
    end
end
predictions = containers.Map();                  % fichiers de resultats attendus, lus une fois

% --- Simulations ---
nouveaux = struct('code', {}, 'nom', {}, 't', {}, 'v', {}, 'consigne', {}, 'mu', {}, 'd', {}, 'iL', {}, ...
                  'K', {}, 'grandeurs', {}, 'duree_calcul_s', {}, 'us_par_pas', {}, 'modele', {});
sources = {};                                    % modele source de chaque entree
banc_iae = [];                                   % IAE du banc (30 ms a la fin), NaN si absent
for i = 1:numel(liste)
    L = liste(i);
    MDL = L.modele;
    sc = charger_scenario(L.code);
    fprintf('\n=== %s : %s, %s (%.0f ms) ===\n', MDL, L.methode, sc.nom, sc.duree * 1e3);
    if bdIsLoaded(MDL)
        if strcmp(get_param(MDL, 'Dirty'), 'on')
            error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', MDL);
        end
        close_system(MDL, 0);
    end
    load_system(fullfile(pwd, [MDL '.slx']));
    try
        % Le modele doit etre celui de Construire_Modeles_Comparaison.m
        stop_attendu = (numel(sc.t) - 0.5) * TC;
        stop = str2double(get_param(MDL, 'StopTime'));
        if ~(abs(stop - stop_attendu) <= 1e-12)
            error('%s : StopTime %s au lieu de %.17g. Reconstruire avec Construire_Modeles_Comparaison.m.', ...
                  MDL, get_param(MDL, 'StopTime'), stop_attendu);
        end
        if ~contains(get_param(MDL, 'InitFcn'), sprintf('charger_scenario(''%s'');', L.code))
            error('%s : le rappel InitFcn ne charge pas %s. Reconstruire avec Construire_Modeles_Comparaison.m.', ...
                  MDL, L.code);
        end
        tw = find_system(MDL, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
        variables = get_param(tw, 'VariableName');
        for v_ = {'cmp_vout', 'cmp_iL', 'cmp_d'}
            if ~any(strcmp(variables, v_{1}))
                error(['%s : enregistrement %s absent (modele fait par un autre script ?). ' ...
                       'Reconstruire avec Construire_Modeles_Comparaison.m.'], MDL, v_{1});
            end
        end
        avec_gains = any(strcmp(variables, 'cmp_K'));
        gains_fixes = NaN(1, 3);                     % gains constants (PSO-PID, ZN), pour la figure
        if FIGURES && ~avec_gains
            gains_fixes = lire_gains_pid(MDL);
        end
        if NEUTRALISER_TO_WORKSPACE_DU_MODELE
            for j = 1:numel(tw)
                if ~startsWith(get_param(tw{j}, 'VariableName'), 'cmp_')
                    set_param(tw{j}, 'Commented', 'on');   % en memoire seulement
                end
            end
        end
        t_debut = tic;
        sortie = sim(MDL, 'ReturnWorkspaceOutputs', 'on');
        duree_calcul = toc(t_debut);
    catch ME
        close_system(MDL, 0);
        rethrow(ME);
    end
    close_system(MDL, 0);                                 % SANS enregistrer

    n_sc = numel(sc.t);
    [t, v] = extraire_enregistrement(sortie, 'cmp_vout');
    [t2, iL] = extraire_enregistrement(sortie, 'cmp_iL');
    [t3, d] = extraire_enregistrement(sortie, 'cmp_d');
    if numel(t) >= n_sc && numel(t2) >= n_sc && numel(t3) >= n_sc
        t = t(1:n_sc); t2 = t2(1:n_sc); t3 = t3(1:n_sc);
        v = v(1:n_sc); iL = iL(1:n_sc); d = d(1:n_sc);
    end
    if ~isequal(numel(t), numel(t2), numel(t3), n_sc) || max(abs([t - t2; t - t3; t - sc.t])) > 1e-9
        error('%s : enregistrements non alignes avec le profil (%d instants attendus, %d obtenus).', MDL, n_sc, numel(t));
    end
    K = [];
    if avec_gains
        [tk, K] = extraire_matrice(sortie, 'cmp_K');
        if size(K, 1) >= n_sc
            K = K(1:n_sc, :);
            tk = tk(1:n_sc);
        end
        if size(K, 1) ~= n_sc || max(abs(tk - sc.t)) > 1e-9
            error('%s : enregistrement des gains non aligne avec le profil.', MDL);
        end
    end
    if any(d < 0.01 - 1e-9 | d > 0.99 + 1e-9)
        error(['%s : le signal qui entre dans le PWM sort de [0.01 ; 0.99] (de %.4f a %.4f). ' ...
               'Il manque la saturation physique devant le PWM.'], MDL, min(d), max(d));
    end
    consigne = 100 + sc.dvref;
    g = grandeurs(v, consigne, d, iL, sc.evenements, TC);
    us = duree_calcul / n_sc * 1e6;
    % mu : sortie du bloc regulateur. Dans les cinq modeles, la sortie du PID
    % (saturation [0.01 ; 0.99] comprise) est l'entree du PWM : mu = d.
    nouveaux(end + 1) = struct('code', sc.code, 'nom', sc.nom, 't', t, 'v', v, 'consigne', consigne, 'mu', d, ...
                               'd', d, 'iL', iL, 'K', K, 'grandeurs', g, 'duree_calcul_s', duree_calcul, ...
                               'us_par_pas', us, 'modele', MDL); %#ok<SAGROW>
    sources{end + 1} = L.source; %#ok<SAGROW>

    fprintf('  temps de calcul : %.1f s (%.1f us par periode du regulateur)\n', duree_calcul, us);
    fprintf('  Simulink : %s\n', resume(g));
    for j = 1:numel(g.evenements)
        ev = g.evenements(j);
        fprintf('    evenement a %6.1f ms : IAE %7.2f mV.s, ecart max %6.2f V, %s\n', ev.t_ms, ev.IAE * 1e3, ...
                ev.e_max_V, texte_retour(ev));
    end
    if ~isempty(K)
        fprintf('  gains finaux [P I D] : [%.5g %.5g %.5g]\n', K(end, :));
    end
    % Comparaison au banc (resultats attendus de la methode)
    if ~isKey(predictions, L.predictions)
        pred = [];
        if isfile(fullfile(pwd, L.predictions))
            pred = jsondecode(fileread(fullfile(pwd, L.predictions)));
            if ~isfield(pred, 'format') || pred.format ~= 2
                warning('%s n''a pas le format attendu : pas de comparaison.', L.predictions);
                pred = [];
            end
        end
        predictions(L.predictions) = pred;
    end
    pred = predictions(L.predictions);
    banc_iae(end + 1) = NaN; %#ok<SAGROW>
    if ~isempty(pred) && isfield(pred.essais, L.code)
        gb = pred.essais.(L.code).grandeurs;
        vb = pred.essais.(L.code).v_toutes_les_ms(:);
        vs = v(1:round(1e-3 / TC):end);
        n_cmp = min(numel(vb), numel(vs));
        banc_iae(end) = gb.IAE;
        fprintf('  banc     : %s\n', resume(gb));
        fprintf(['  IAE de 30 ms a la fin : Simulink %.3f, banc %.3f mV.s (%+.1f %%) ; ecart Simulink - banc sur ' ...
                 'Vout (toutes les ms) : %.1f mV efficace, %.1f mV au plus\n'], g.IAE * 1e3, gb.IAE * 1e3, ...
                (g.IAE / gb.IAE - 1) * 100, sqrt(mean((vs(1:n_cmp) - vb(1:n_cmp)).^2)) * 1e3, ...
                max(abs(vs(1:n_cmp) - vb(1:n_cmp))) * 1e3);
        if ~isempty(K) && isfield(pred.essais.(L.code), 'K_toutes_les_ms')
            Kb = pred.essais.(L.code).K_toutes_les_ms;
            Ks = K(1:round(1e-3 / TC):end, :);
            n_k = min(size(Kb, 1), size(Ks, 1));
            fprintf('  gains : ecart relatif maximal Simulink - banc (toutes les ms) %.2g\n', ...
                    max(max(abs(Ks(1:n_k, :) - Kb(1:n_k, :)) ./ abs(Kb(1:n_k, :)))));
        end
        evb = gb.evenements;
        if iscell(evb)
            evb = [evb{:}];
        end
        for j = 1:min(numel(evb), numel(g.evenements))
            fprintf('    evenement a %6.1f ms : IAE %7.2f / %7.2f mV.s ; retour : %s / %s (Simulink / banc)\n', ...
                    g.evenements(j).t_ms, g.evenements(j).IAE * 1e3, evb(j).IAE * 1e3, ...
                    texte_retour(g.evenements(j)), texte_retour(evb(j)));
        end
    end
    % Figure (apres la mesure du temps ; n'entre dans aucune grandeur)
    if FIGURES
        try
            fig = tracer_figure(MDL, L.methode, sc, t, v, consigne, d, iL, K, gains_fixes, K_ZN, MARGE_ORDONNEES);
            fichiers_fig = enregistrer_figure(fig, fullfile(pwd, DOSSIER_FIGURES), MDL);
            fprintf('  figure : %s\n', strjoin(fichiers_fig, ' et '));
            if FERMER_FIGURES
                close(fig);
            end
        catch ME
            warning('%s : figure non produite (%s). Les resultats ne sont pas touches.', MDL, ME.message);
        end
    end
end

% --- Fichiers : un par modele source, au format de Simuler_<Methode>.m ---
for s = unique(sources, 'stable')
    MODELE = s{1}; %#ok<NASGU>
    fichier = fullfile(pwd, ['resultats_' s{1} '_comparaison.mat']);
    neufs = nouveaux(strcmp(sources, s{1}));
    resultats = neufs;
    if isfile(fichier)                                    % garder les scenarios non resimules
        try
            ancien = load(fichier, 'resultats');
            garde = ancien.resultats(~ismember({ancien.resultats.code}, {neufs.code}));
            resultats = [garde(:)', neufs(:)'];
            rang = zeros(1, numel(resultats));            % ordre de CODES, les autres a la suite
            for q = 1:numel(resultats)
                r_ = find(strcmp(CODES, resultats(q).code), 1);
                if isempty(r_)
                    r_ = numel(CODES) + q;
                end
                rang(q) = r_;
            end
            [~, o] = sort(rang);
            resultats = resultats(o);
            if ~isempty(garde)
                fprintf('\n%s : scenarios gardes de la simulation precedente : %s.\n', fichier, ...
                        strjoin({garde.code}, ', '));
            end
        catch ME
            warning('%s illisible ou d''un autre format (%s) : remplace.', fichier, ME.message);
            resultats = neufs;
        end
    end
    save(fichier, 'resultats', 'MODELE', '-v7.3');
    fprintf('\n%s ecrit (%s).\n', fichier, strjoin({resultats.code}, ', '));
end

% --- Tableau recapitulatif ---
fprintf('\n==========================================================================================\n');
fprintf(' Recapitulatif (Simulink). Les metriques completes : COMPARAISON/metriques/Metriques_Simulink.m,\n');
fprintf(' reglage RESULTATS = ''comparaison''.\n');
fprintf('==========================================================================================\n');
fprintf('  %-14s %-16s %10s %10s %14s %14s %8s\n', 'modele', 'methode', 'calcul s', 'us / Tc', ...
        'IAE30-fin mV.s', 'banc mV.s', 'ecart %');
for i = 1:numel(nouveaux)
    g = nouveaux(i).grandeurs;
    fprintf('  %-14s %-16s %10.1f %10.1f %14.3f %14.3f %8.2f\n', nouveaux(i).modele, liste(i).methode, ...
            nouveaux(i).duree_calcul_s, nouveaux(i).us_par_pas, g.IAE * 1e3, banc_iae(i) * 1e3, ...
            (g.IAE / banc_iae(i) - 1) * 100);
end


%% ===================== Fonctions =====================================

function o = grandeurs(v, consigne, d, iL, evenements, Te)
    % Memes definitions que la fonction grandeurs de banc_commun.py (copie
    % sans changement de Simuler_ELM_PID.m). Les numeros de pas commencent a
    % 0 (comme en Python) ; l'indice MATLAB du pas n est n + 1.
    N = numel(v);
    n = (0:N-1)';
    e = consigne - v;
    k30 = round(0.03 / Te);
    evts = round(evenements(:)' / Te);
    k_dem = min([round(0.05 / Te), evts]);                % fin du demarrage (exclue)
    % Demarrage
    k10 = find(v >= 10, 1);
    k90 = find(v >= 90, 1);
    if isempty(k10) || isempty(k90)
        o.t_montee_ms = NaN;
    else
        o.t_montee_ms = (k90 - k10) * Te * 1e3;
    end
    o.depassement_pct = max(0, (max(v(1:k_dem)) - 100) / 100 * 100);
    hors = find(abs(e(1:k_dem)) > 1);
    if isempty(hors)
        o.t_etab_ms = 0;
    else
        o.t_etab_ms = hors(end) * Te * 1e3;               % (indice Python + 1) x Te
    end
    o.erreur_moy_dem_V = mean(e(k30 + 1:k_dem));
    o.ondulation_V = max(v(k30 + 1:k_dem)) - min(v(k30 + 1:k_dem));
    o.iL_max_dem_A = max(iL(1:k_dem));                    % courant crete du demarrage (observe)
    o.iL_ondulation_A = max(iL(k30 + 1:k_dem)) - min(iL(k30 + 1:k_dem));   % ondulation de courant en regime
    % De 30 ms a la fin
    w = n >= k30;
    o.IAE = sum(abs(e(w))) * Te;
    o.ISE = sum(e(w).^2) * Te;
    o.ITAE = sum(n(w) * Te .* abs(e(w))) * Te;
    o.e_max_V = max(abs(e(w)));
    o.e_eff_V = sqrt(mean(e(w).^2));
    % Evenements
    o.evenements = struct('t_ms', {}, 'IAE', {}, 'e_max_V', {}, 't_retour_ms', {}, 'reste_dans_bande', {}, ...
                          'revenu', {}, 'e_crete_a_crete_fin_V', {}, 'iL_max_A', {}, 'iL_min_A', {});
    bornes = [evts, N];
    n_fin = round(0.010 / Te);                            % 10 dernieres ms de chaque fenetre (2200 instants)
    for j = 1:numel(evts)
        fen = n >= evts(j) & n < bornes(j + 1);           % de l'evenement au suivant
        fin = n >= max(evts(j), bornes(j + 1) - n_fin) & n < bornes(j + 1);   % ses 10 dernieres ms
        hors = find(fen & abs(e) > 1);                    % indices MATLAB des pas hors bande
        revenu = max(abs(e(fin))) <= 1;                   % dans la bande sur les 10 dernieres ms
        if isempty(hors)
            t_ret = 0;                                    % jamais sorti de la bande
        elseif revenu
            t_ret = (hors(end) - evts(j)) * Te * 1e3;     % (indice Python - kj + 1) x Te
        else
            t_ret = NaN;                                  % pas revenu avant l'evenement suivant
        end
        o.evenements(j) = struct('t_ms', evts(j) * Te * 1e3, 'IAE', sum(abs(e(fen))) * Te, ...
                                 'e_max_V', max(abs(e(fen))), 't_retour_ms', t_ret, ...
                                 'reste_dans_bande', isempty(hors), 'revenu', revenu, ...
                                 'e_crete_a_crete_fin_V', max(e(fin)) - min(e(fin)), ...
                                 'iL_max_A', max(iL(fen)), 'iL_min_A', min(iL(fen)));
    end
    % Activite de la commande et conduction discontinue
    dw = d(w);
    o.dd_moyen = mean(abs(diff(dw)));
    o.butee_pct = mean(dw <= 0.01 + 1e-12 | dw >= 0.99 - 1e-12) * 100;
    o.d_ecart_type = std(dw, 1);                          % ecart type "population", comme numpy
    o.dcm_pct = mean(iL(w) < 0.01) * 100;
end

function s = resume(g)
    % Une ligne lisible des grandeurs principales (memes champs que le banc).
    s = sprintf(['montee %.2f ms, depassement %.2f %%, etabli a %.2f ms ; IAE %.2f mV.s, ecart max %.2f V, ' ...
                 'efficace %.1f mV ; |dd| %.4f, butee %.1f %%, DCM %.1f %%'], g.t_montee_ms, g.depassement_pct, ...
                g.t_etab_ms, g.IAE * 1e3, g.e_max_V, g.e_eff_V * 1e3, g.dd_moyen, g.butee_pct, g.dcm_pct);
    if isfield(g, 'iL_max_dem_A')
        s = sprintf('%s ; iL crete au demarrage %.1f A', s, g.iL_max_dem_A);
    end
end

function s = texte_retour(ev)
    % Etat de retour dans la bande d'un evenement (Simulink ou banc).
    if ev.reste_dans_bande
        s = 'reste dans la bande';
    elseif ev.revenu
        s = sprintf('retour en %.2f ms', ev.t_retour_ms);
    else
        s = 'PAS REVENU';
    end
end

function [t, M] = extraire_matrice(sortie, nom)
    % Lit un enregistrement To Workspace d'un signal vectoriel : une ligne
    % par instant, une colonne par composante.
    if ~any(strcmp(sortie.who, nom))
        error('Enregistrement "%s" absent des resultats de simulation.', nom);
    end
    donnee = sortie.get(nom);
    if isa(donnee, 'timeseries')
        t = donnee.Time;
        D = donnee.Data;
    else
        t = donnee.time;
        D = donnee.signals.values;
    end
    t = t(:);
    D = double(squeeze(D));
    if size(D, 1) ~= numel(t)                             % composantes en lignes : on transpose
        D = D.';
    end
    M = D;
end

function [t, v] = extraire_enregistrement(sortie, nom)
    % Lit un enregistrement To Workspace dans l'objet renvoye par sim().
    if ~any(strcmp(sortie.who, nom))
        error('Enregistrement "%s" absent des resultats de simulation.', nom);
    end
    donnee = sortie.get(nom);
    if isa(donnee, 'timeseries')
        t = donnee.Time;
        v = squeeze(donnee.Data);
    else
        t = donnee.time;
        v = squeeze(donnee.signals.values);
    end
    t = t(:);
    v = double(v(:));
end

function fig = tracer_figure(mdl, methode, sc, t, v, consigne, d, iL, K, gains_fixes, K_ZN, marge)
    % Une figure par modele :
    % Vout et consigne, rapport cyclique, iL, gains / ZN (gains constants :
    % droites portant leur valeur et, sous la tuile, le texte des gains ;
    % texte seul si les gains ne se lisent pas). Abscisses en ms sur toute
    % la duree ; ordonnees de chaque tuile fixees par fixer_ordonnees (marge
    % de chaque cote).
    t_ms = t(:) * 1e3;
    evts_ms = reshape(sc.evenements, 1, []) * 1e3;
    avec_courbes_gains = ~isempty(K);
    gains_lus = ~avec_courbes_gains && all(isfinite(gains_fixes));   % gains constants lus dans le bloc
    position = [60 40 1100 820];
    if gains_lus
        % Deux lignes de texte sous la tuile 4 : figure plus haute, posee
        % au-dessus d'une barre des taches et sans depasser le haut de l'ecran.
        ecran = get(groot, 'ScreenSize');
        position = [60 70 1100 max(650, min(900, ecran(4) - 140))];
    end
    fig = figure('Name', mdl, 'Color', 'w', 'Position', position);
    tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s : %s, scenario %s (%s)', mdl, methode, sc.code, sc.nom), 'Interpreter', 'none');
    ax = gobjects(1, 4);
    for p = 1:4
        ax(p) = nexttile(tl);
        if p <= 3 || avec_courbes_gains || gains_lus
            hold(ax(p), 'on');
            grid(ax(p), 'on');
            box(ax(p), 'on');
            for k = 1:numel(evts_ms)
                xline(ax(p), evts_ms(k), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
    end
    % 1. Vout et consigne
    plot(ax(1), t_ms, v, 'LineWidth', 0.8, 'DisplayName', 'Vout');
    plot(ax(1), t_ms, consigne, 'k--', 'LineWidth', 0.8, 'DisplayName', 'consigne');
    legend(ax(1), 'Location', 'northoutside', 'Orientation', 'horizontal');
    ylabel(ax(1), 'Vout (V)');
    fixer_ordonnees(ax(1), [v(:); consigne(:)], marge);
    % 2. Rapport cyclique
    plot(ax(2), t_ms, d, 'LineWidth', 0.6);
    ylabel(ax(2), 'rapport cyclique');
    fixer_ordonnees(ax(2), d, marge);
    % 3. Courant de la bobine
    plot(ax(3), t_ms, iL, 'LineWidth', 0.6);
    ylabel(ax(3), 'iL (A)');
    fixer_ordonnees(ax(3), iL, marge);
    % 4. Gains rapportes a ceux de Ziegler-Nichols
    if avec_courbes_gains
        Kr = K ./ reshape(K_ZN, 1, []);
        plot(ax(4), t_ms, Kr(:, 1), '-', 'LineWidth', 1.0, 'DisplayName', 'P');
        plot(ax(4), t_ms, Kr(:, 2), '--', 'LineWidth', 1.0, 'DisplayName', 'I');
        plot(ax(4), t_ms, Kr(:, 3), ':', 'LineWidth', 1.0, 'DisplayName', 'D');
        legend(ax(4), 'Location', 'northoutside', 'Orientation', 'horizontal');
        ylabel(ax(4), 'gains / ZN');
        fixer_ordonnees(ax(4), Kr(:), marge);
        xlabel(ax(4), 'temps (ms)');
        lies = ax;
    elseif gains_lus
        % Gains constants : une droite
        % pleine par valeur distincte des rapports P, I, D a ZN, sans
        % marqueur, couleurs par defaut dans l'ordre P, I, D comme pour les
        % methodes adaptatives. La valeur est ecrite sur chaque droite (gain
        % lu et rapport a ZN, memes formats que texte_gains_fixes). Rapports
        % egaux a 5e-4 pres en relatif, la moitie du dernier chiffre du
        % rapport affiche en %.3f (ZN : tous valent 1, meme si les gains du
        % bloc sont arrondis) : une seule droite et une seule etiquette pour
        % le groupe. Etiquettes a 15, 45 et 75 % de la duree (de la plus
        % basse a la plus haute droite) :
        % au-dessus de la droite la plus haute, au-dessous de la plus basse,
        % et pour celle du milieu du cote oppose a sa plus proche voisine ;
        % deux droites proches n'ont donc jamais leurs etiquettes du meme
        % cote. Ordonnees : marge de 35 % de chaque cote, elargie si besoin
        % jusqu'a ce que les etiquettes tiennent dans la tuile sans toucher
        % les bornes. Graduations aux valeurs des droites si elles sont
        % assez espacees (10 % de l'echelle), sinon graduations ordinaires.
        % Le texte des gains reste sous la tuile.
        g = reshape(gains_fixes, 1, []);
        Kr = g ./ reshape(K_ZN, 1, []);
        noms = {'P', 'I', 'D'};
        couleurs = get(ax(4), 'ColorOrder');
        groupe = zeros(1, 3);                             % groupes de rapports egaux
        n_gr = 0;
        for q = 1:3
            if groupe(q) == 0
                n_gr = n_gr + 1;
                egaux = groupe == 0 & abs(Kr - Kr(q)) <= 5e-4 * max(abs(Kr), abs(Kr(q)));
                groupe(egaux) = n_gr;
            end
        end
        val = zeros(1, n_gr);
        etiq = cell(1, n_gr);
        coul = zeros(n_gr, 3);
        for k = 1:n_gr
            m = find(groupe == k);
            val(k) = Kr(m(1));
            coul(k, :) = couleurs(m(1), :);
            if isscalar(m)
                etiq{k} = sprintf('%s = %.5g (%.3f x ZN)', noms{m}, g(m), val(k));
            else
                detail = arrayfun(@(j) sprintf('%s = %.5g', noms{j}, g(j)), m, 'UniformOutput', false);
                etiq{k} = sprintf('%s = %.3f x ZN (%s)', strjoin(noms(m), ', '), val(k), strjoin(detail, ', '));
            end
            plot(ax(4), [t_ms(1), t_ms(end)], [val(k), val(k)], '-', 'LineWidth', 2, ...
                 'Color', coul(k, :), 'DisplayName', strjoin(noms(m), ' = '));
        end
        legend(ax(4), 'Location', 'northoutside', 'Orientation', 'horizontal');
        ylabel(ax(4), 'gains / ZN');
        xlabel(ax(4), [{'temps (ms)'}, texte_gains_fixes(gains_fixes, K_ZN)], 'Interpreter', 'none');
        fixer_ordonnees(ax(4), val, max(marge, 0.35));
        % Cote et abscisse de chaque etiquette
        [vs, ordre] = sort(val);
        dessus = true(1, n_gr);
        frac_x = [0.15, 0.45, 0.75];                      % part de la duree
        x_et = zeros(1, n_gr);
        for r = 1:n_gr
            k = ordre(r);
            x_et(k) = t_ms(1) + frac_x(r) * (t_ms(end) - t_ms(1));
            if r == 1 && n_gr > 1
                dessus(k) = false;                        % plus basse : au-dessous
            elseif r > 1 && r < n_gr                      % milieu : loin de la plus proche voisine
                dessus(k) = vs(r + 1) - vs(r) >= vs(r) - vs(r - 1);
            end
        end
        h_et = gobjects(1, n_gr);
        for k = 1:n_gr
            if dessus(k)
                alignement = 'bottom';
            else
                alignement = 'top';
            end
            h_et(k) = text(ax(4), x_et(k), val(k), etiq{k}, 'Color', coul(k, :), 'FontSize', 9, ...
                           'HorizontalAlignment', 'left', 'VerticalAlignment', alignement, ...
                           'BackgroundColor', 'w', 'Margin', 1, 'Interpreter', 'none');
        end
        % Ecart droite-etiquette de 4 % de l'echelle ; ordonnees elargies
        % tant qu'une etiquette approche une borne a moins de 3 %.
        for iter = 1:5
            yl = ylim(ax(4));
            for k = 1:n_gr
                set(h_et(k), 'Position', [x_et(k), val(k) + (2 * dessus(k) - 1) * 0.04 * diff(yl), 0]);
            end
            drawnow;
            ext = zeros(n_gr, 4);
            for k = 1:n_gr
                ext(k, :) = get(h_et(k), 'Extent');
            end
            bas = min(ext(:, 2));
            haut = max(ext(:, 2) + ext(:, 4));
            garde = 0.03 * diff(yl);
            if bas >= yl(1) + garde && haut <= yl(2) - garde
                break;
            end
            ylim(ax(4), [min(yl(1), bas - 2 * garde), max(yl(2), haut + 2 * garde)]);
        end
        if n_gr == 1 || min(diff(vs)) >= 0.1 * diff(ylim(ax(4)))
            set(ax(4), 'YTick', vs, 'YTickLabel', arrayfun(@(y) sprintf('%.3f', y), vs, 'UniformOutput', false));
        end
        lies = ax;
    else
        axis(ax(4), 'off');
        text(ax(4), 0.5, 0.5, texte_gains_fixes(gains_fixes, K_ZN), 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11, ...
             'Interpreter', 'none');
        xlabel(ax(3), 'temps (ms)');
        lies = ax(1:3);
    end
    linkaxes(lies, 'x');
    xlim(ax(1), [t_ms(1), t_ms(end)]);                    % toute la duree de l'essai
end

function fixer_ordonnees(ax, y, marge)
    % Ordonnees de la tuile : [min - m ; max + m] des donnees tracees, avec
    % m = marge x etendue. Courbe plate (etendue nulle) : etendue prise a
    % 10 % de la valeur (1 si la valeur est nulle), pour que la courbe soit
    % au milieu de la tuile et loin des bornes.
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

function txt = texte_gains_fixes(g, K_ZN)
    % Texte de la tuile 4 quand la methode n'a pas de gains variables (deux
    % lignes, en cellule).
    if all(isfinite(g))
        txt = {sprintf('Gains constants (pas d''adaptation) : P = %.5g, I = %.5g, D = %.5g', g), ...
               sprintf('soit %.3f / %.3f / %.3f fois ceux de Ziegler-Nichols', g ./ reshape(K_ZN, 1, []))};
    else
        txt = {'Gains constants (pas d''adaptation) : pas de courbe de gains.', ...
               '(valeurs non lues dans le bloc "PID Controller")'};
    end
end

function g = lire_gains_pid(mdl)
    % Gains P, I, D du bloc "PID Controller" (modeles a gains constants :
    % PSO-PID, Ziegler-Nichols). NaN si le bloc ou une valeur ne se lit pas ;
    % sans consequence sur la simulation.
    g = NaN(1, 3);
    bloc = [mdl '/PID Controller'];
    noms = {'P', 'I', 'D'};
    try
        if getSimulinkBlockHandle(bloc) == -1
            return;
        end
        for q = 1:3
            s = get_param(bloc, noms{q});
            x = str2double(s);
            if isnan(x)
                try
                    x = evalin('base', s);               % gain donne par une variable
                catch
                    x = NaN;
                end
            end
            if isnumeric(x) && isscalar(x)
                g(q) = double(x);
            end
        end
    catch
        g = NaN(1, 3);
    end
end

function fichiers = enregistrer_figure(fig, dossier, nom)
    % <nom>.png (200 dpi) et <nom>.fig dans le dossier (cree s'il manque).
    if ~isfolder(dossier)
        mkdir(dossier);
    end
    fichiers = {fullfile(dossier, [nom '.png']), fullfile(dossier, [nom '.fig'])};
    exportgraphics(fig, fichiers{1}, 'Resolution', 200);
    savefig(fig, fichiers{2});
end
