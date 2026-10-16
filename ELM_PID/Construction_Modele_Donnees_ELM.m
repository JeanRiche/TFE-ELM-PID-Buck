%% CONSTRUCTION_MODELE_DONNEES_ELM.m
%
% OBJECTIF
% --------
% Construire Donnees_ELM.slx, le modele qui enregistre le convertisseur en
% boucle ouverte pour l'apprentissage de l'ELM. C'est une copie de
% Buck_Commun.slx dans laquelle le bloc
% "PID Controller" est remplace par une source de rapport cyclique impose.
%
% CE QUE LE SCRIPT CHANGE, ET RIEN D'AUTRE
% -----------------------------------------
%   1. Le bloc "PID Controller" est supprime. A sa place, "Commande
%      Excitation" (From Workspace) lit le profil sc_dexc a la periode du
%      regulateur Tc = 1/(22000*10) s, sans interpolation, et alimente
%      l'entree du PWM comme le faisait le PID.
%   2. La sortie de Sum1 (l'erreur), qui n'alimente plus rien, va dans
%      "Terminaison Erreur" (Terminator), pour qu'aucune sortie ne reste en
%      l'air.
%   3. Le rappel InitFcn charge aussi un rapport cyclique de 0.5 si aucun
%      profil sc_dexc n'existe (le modele reste utilisable avec Run).
% Le circuit de puissance, les sources commandees (Vin, charge
% electronique), la mesure, le PWM, le powergui et la configuration ne
% sont pas modifies. Le rapport cyclique impose reste dans [0.05 ; 0.95] :
% la saturation du PID n'a pas a etre remplacee.
%
% LES AUTO-TESTS (apres la sauvegarde)
% -------------------------------------
%   T1  demarrage de ELM_A1 (40 ms a d = 0.5, 5 ohms, 200 V) : tension
%       moyenne par fenetre de 0.5 ms comparee au banc commun ;
%   T2  les 100 premieres ms de ELM_A1 (premiers paliers de d, rampes de
%       Vin et de charge) : meme comparaison ;
%       ecart maximal tolere 0.05 V (le banc a predit Simulink a quelques
%       mV pres sur les essais communs) ;
%   T3  le rapport cyclique qui entre dans le PWM est exactement sc_dexc.
%
% ORDRE D'EXECUTION (donnees de l'ELM-PID)
% -----------------------------------------
%   1. excitation_donnees_elm.py ; 2. ce script ; 3. Generer_Donnees_ELM.m.
%
% Prerequis dans le dossier courant MATLAB : Buck_Commun.slx,
% charger_scenario.m, scenario_S1.mat, scenario_ELM_A1.mat et
% predictions_banc_donnees_elm.mat.
% Duree : une a trois minutes (auto-tests compris). Compatible R2024a.
% ---------------------------------------------------------------------

clear; clc;

SOURCE = 'Buck_Commun';                           % modele commun (jamais modifie)
MDL    = 'Donnees_ELM';                           % modele construit
TC_TXT = '1/(22000*10)';                          % periode du regulateur (texte exact)
TC     = 1/(22000*10);                            % la meme, en nombre (s)
NF     = 110;                                     % instants Tc par fenetre de 0.5 ms

% --- Prerequis ---
for f = {[SOURCE '.slx'], 'charger_scenario.m', 'scenario_S1.mat', 'scenario_ELM_A1.mat', ...
         'predictions_banc_donnees_elm.mat'}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
end
if numel(which('charger_scenario', '-all')) ~= 1
    error('Il doit exister exactement un charger_scenario.m sur le chemin MATLAB.');
end
if bdIsLoaded(SOURCE) && strcmp(get_param(SOURCE, 'Dirty'), 'on')
    error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', SOURCE);
end
for m = {MDL, SOURCE}
    if bdIsLoaded(m{1})
        close_system(m{1}, 0);                    % fermer sans enregistrer
    end
end
charger_scenario('S1');                           % variables sc_* presentes avant le chargement
assignin('base', 'sc_dexc', [0 0.5]);

%% 1. Copie fraiche
chemin = fullfile(pwd, [MDL '.slx']);
if isfile(chemin)
    delete(chemin);
    fprintf('Ancien %s.slx supprime (on repart toujours d''une copie fraiche).\n', MDL);
end
copyfile(fullfile(pwd, [SOURCE '.slx']), chemin);
if numel(which([MDL '.slx'], '-all')) ~= 1
    error('Il doit exister exactement un %s.slx sur le chemin MATLAB.', MDL);
end
load_system(chemin);
fprintf('Copie de %s.slx chargee sous le nom %s.\n', SOURCE, MDL);

%% 2. Remplacement du PID par la commande imposee
pid = [MDL '/PID Controller'];
if getSimulinkBlockHandle(pid) == -1
    close_system(MDL, 0);
    error('Bloc introuvable : %s', pid);
end
lh = get_param(pid, 'LineHandles');
[src_err, ~] = extremites_ligne(lh.Inport(1));    % ce qui alimente le PID (Sum1, sortie 1)
[~, dst_cmd] = extremites_ligne(lh.Outport(1));   % ce que le PID alimente (entree du PWM)
if ~strcmp(src_err.bloc, [MDL '/Sum1']) || numel(dst_cmd) ~= 1 || ...
        ~strcmp(dst_cmd(1).bloc, [MDL '/' sprintf('PWM Generator\n(DC-DC)')])
    close_system(MDL, 0);
    error('Le PID devrait etre alimente par Sum1 et alimenter seulement le PWM.');
end
pos = get_param(pid, 'Position');
delete_line(lh.Inport(1));                        % lignes du PID seulement
delete_line(lh.Outport(1));
delete_block(pid);

cmd = [MDL '/Commande Excitation'];
add_block('simulink/Sources/From Workspace', cmd, 'Position', pos);
set_param(cmd, 'VariableName', 'sc_dexc', 'SampleTime', TC_TXT, 'Interpolate', 'off', ...
          'OutputAfterFinalValue', 'Holding final value');
relier(MDL, cmd, 1, dst_cmd(1).bloc, dst_cmd(1).port);

term = [MDL '/Terminaison Erreur'];
add_block('simulink/Sinks/Terminator', term, 'Position', [pos(1), pos(4) + 30, pos(1) + 20, pos(4) + 50]);
relier(MDL, src_err.bloc, src_err.port, term, 1);
fprintf('PID remplace par "Commande Excitation" (sc_dexc, %s s) ; Sum1 vers "Terminaison Erreur".\n', TC_TXT);

%% 3. Rappel InitFcn, hygiene et sauvegarde
set_param(MDL, 'InitFcn', ['if ~exist(''sc_R0'', ''var''), charger_scenario(''S1''); end; ' ...
                           'if ~exist(''sc_dexc'', ''var''), sc_dexc = [0 0.5]; end']);
set_param(MDL, 'ReturnWorkspaceOutputs', 'off');
save_system(MDL);
close_system(MDL, 0);
fprintf('%s.slx sauvegarde.\n', MDL);

%% 4. Auto-tests
fprintf('\n=== Auto-tests ===\n');
P = load(fullfile(pwd, 'predictions_banc_donnees_elm.mat'));
v_banc = P.v_fen_ELM_A1(:);                       % tension moyenne par fenetre, banc commun
for essai = {{'T1 (demarrage, 40 ms)', 0.040}, {'T2 (100 ms de ELM_A1)', 0.100}}
    nom = essai{1}{1};
    duree = essai{1}{2};
    sc = charger_scenario('ELM_A1');
    dexc = charger_commande('ELM_A1', sc);
    r = simuler_court(MDL, duree);
    n_fen = floor(numel(r.v) / NF);
    v_fen = mean(reshape(r.v(1:n_fen * NF), NF, n_fen), 1)';
    ecart = max(abs(v_fen - v_banc(1:n_fen)));
    fprintf('  %s : %d fenetres, ecart maximal au banc %.1f mV (tension de %.1f a %.1f V)\n', ...
            nom, n_fen, ecart * 1e3, min(v_fen), max(v_fen));
    if ~(ecart <= 0.05)
        error('%s : ecart au banc de %.3f V, au-dela de 0.05 V.', nom, ecart);
    end
    ecart_d = max(abs(r.d - dexc(1:numel(r.d))));
    if ~(ecart_d <= 1e-12)
        error('T3 : le rapport cyclique applique differe de sc_dexc de %.3g.', ecart_d);
    end
end
fprintf('  T3 reussi : le rapport cyclique applique est exactement sc_dexc.\n');
fprintf(['\n%s.slx est construit et ses auto-tests sont passes.\n' ...
         'Lancer ensuite Generer_Donnees_ELM.m.\n'], MDL);


%% ===================== Fonctions =====================================

function dexc = charger_commande(code, sc)
    % Ecrit sc_dexc dans le workspace de base avec les memes instants
    % decales d'un demi-pas que charger_scenario (voir son en-tete), et
    % renvoie le profil brut (un rapport cyclique par instant k*Tc).
    m = load(fullfile(pwd, ['scenario_' code '.mat']), 'dexc');
    dexc = double(m.dexc(:));
    n = numel(sc.t);
    if numel(dexc) ~= n
        error('scenario_%s.mat : %d valeurs de dexc pour %d instants.', code, numel(dexc), n);
    end
    td = [0; ((1:n-1)' - 0.5) * sc.Te];
    assignin('base', 'sc_dexc', [td, dexc]);
end

function r = simuler_court(mdl, duree)
    % Simulation avec enregistrements temporaires (Vout, rapport cyclique
    % applique au PWM), modele ferme SANS enregistrer, meme en cas d'erreur.
    load_system(fullfile(pwd, [mdl '.slx']));
    ajouter_enregistrement(mdl, [mdl '/Log Vout tmp'], 'tmp_v', [mdl '/Vout'], 1);
    ajouter_enregistrement(mdl, [mdl '/Log d tmp'], 'tmp_d', [mdl '/Commande Excitation'], 1);
    blocs = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
    for i = 1:numel(blocs)                        % journaux d'origine neutralises en memoire
        if ~startsWith(get_param(blocs{i}, 'VariableName'), 'tmp_')
            set_param(blocs{i}, 'Commented', 'on');
        end
    end
    try
        sortie = sim(mdl, 'StopTime', sprintf('%.17g', duree - 0.5 * 1/(22000*10)), ...
                     'ReturnWorkspaceOutputs', 'on', 'SaveTime', 'off');
    catch ME
        close_system(mdl, 0);
        rethrow(ME);
    end
    close_system(mdl, 0);
    [~, r.v] = extraire_enregistrement(sortie, 'tmp_v');
    [~, r.d] = extraire_enregistrement(sortie, 'tmp_d');
end

function [src, dst] = extremites_ligne(ligne)
    % Source et destinations d'une ligne de signal (les -1
    % sont filtres).
    h_sb = get_param(ligne, 'SrcBlockHandle');
    h_sp = get_param(ligne, 'SrcPortHandle');
    src = struct('bloc', getfullname(h_sb), 'port', get_param(h_sp, 'PortNumber'));
    h_db = get_param(ligne, 'DstBlockHandle');
    h_dp = get_param(ligne, 'DstPortHandle');
    dst = struct('bloc', {}, 'port', {});
    for i = 1:numel(h_db)
        if h_db(i) ~= -1 && h_dp(i) ~= -1
            dst(end+1) = struct('bloc', getfullname(h_db(i)), 'port', get_param(h_dp(i), 'PortNumber')); %#ok<AGROW>
        end
    end
end

function relier(mdl, bloc_src, port_src, bloc_dst, port_dst)
    % Ligne de signal entre deux ports, handles demandes au moment meme.
    ph_src = get_param(bloc_src, 'PortHandles');
    ph_dst = get_param(bloc_dst, 'PortHandles');
    add_line(mdl, ph_src.Outport(port_src), ph_dst.Inport(port_dst), 'autorouting', 'on');
end

function ajouter_enregistrement(mdl, chemin, variable, bloc_src, port_src)
    % Bloc To Workspace temporaire a Tc (texte exact), MaxDataPoints = inf.
    pos = get_param(bloc_src, 'Position');
    add_block('simulink/Sinks/To Workspace', chemin, 'Position', [pos(3)+120, pos(2)+40, pos(3)+200, pos(2)+65]);
    set_param(chemin, 'VariableName', variable, 'SaveFormat', 'Timeseries', ...
              'SampleTime', '1/(22000*10)', 'MaxDataPoints', 'inf');
    relier(mdl, bloc_src, port_src, chemin, 1);
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
