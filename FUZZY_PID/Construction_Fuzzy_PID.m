%% CONSTRUCTION_FUZZY_PID.m
%
% OBJECTIF
% --------
% Construire Buck_Commun_Fuzzy_PID.slx, le modele du Fuzzy-PID de Zhao,
% Tomizuka et Isaka (1993) : une copie de Buck_Commun.slx ou le bloc "PID Controller" de Ziegler-Nichols est GARDE, mais
% recoit P, I et D de l'ordonnanceur flou au lieu de ses champs internes.
%
% CE QUE LE SCRIPT CHANGE, ET RIEN D'AUTRE
% -----------------------------------------
%   1. Bloc "PID Controller" : seul le reglage "Source" des parametres passe
%      de internal a external (ControllerParametersSource). Tout le reste du
%      bloc est inchange : forme Parallel, Forward Euler, saturation
%      [0.01 ; 0.99], clamping, conditions initiales 0.5 et 0.01. Le bloc
%      gagne alors des entrees P, I, D et N, reperees par leur nom (les
%      Inport sous le masque), pas par leur numero. La liaison Sum1 -> PID
%      est refaite sur l'entree "u", ou qu'elle soit.
%   2. "Echantillonneur Erreur Fuzzy" (Zero-Order Hold a 1/(22000*10) s),
%      alimente par Sum1 (l'erreur, comme le PID).
%   3. "Ordonnanceur Flou Zhao" (MATLAB System, classe
%      ordonnanceur_flou_zhao, execution interpretee) : erreur -> [P I D].
%   4. "Demux Gains Fuzzy" (Demux a 3 sorties) -> entrees P, I, D du PID.
%   5. "Constante N Filtre" (Constant) -> entree N du PID, valeur ecrite
%      comme dans le champ N du bloc d'origine (64122.9).
% Le circuit de puissance, les sources commandees, la mesure, la consigne,
% le PWM, le powergui, la configuration et le rappel InitFcn ne changent
% pas. La sortie du PID va toujours au PWM.
%
% LES AUTO-TESTS (apres la sauvegarde ; modele ferme sans enregistrer)
% --------------------------------------------------------------------
%   T0  20 ms de S1, ordonnanceur force aux gains de Ziegler-Nichols
%       (GainsFixes) : Vout toutes les ms comparee au PID classique du banc
%       (predictions_banc_pid_classique.json), 0.1 V au plus. Controle le
%       branchement des entrees externes (ordre P, I, D, N ; entree u).
%   T1  20 ms de S1, Fuzzy-PID : Vout comparee au banc
%       (predictions_banc_fuzzy_pid.json) : 0.5 V au plus sur les 20 ms, et
%       0.01 V au plus de 10 a 20 ms ; gains affiches.
%       Le regulateur Simulink calcule exactement ce que suppose le banc
%       (3.5e-15 sur la commande) ; l'ecart de 1 mV laisse par le circuit au
%       debut est amplifie a 0.2 V entre 2 et 6 ms (apres le pic, quand les
%       gains bougent vite), puis s'eteint (0.7 mV apres 10 ms).
%   T2  20 ms de S8b (98 ohms) : meme comparaison, 0.5 V au plus.
% Si T0 passe et T1 echoue : le bloc PID ne traite pas des gains variables
% comme le banc le suppose.
%
% ORDRE D'EXECUTION
% -----------------
%   1. banc_fuzzy_pid.py ; 2. Tester_Fuzzy_PID_Rejeu.m ; 3. ce script ;
%   4. Simuler_Fuzzy_PID.m.
%
% Prerequis dans le dossier courant MATLAB : Buck_Commun.slx, charger_scenario.m, scenario_S1.mat, scenario_S8b.mat,
% ordonnanceur_flou_zhao.m, fuzzy_pid_reglages.mat,
% predictions_banc_fuzzy_pid.json et predictions_banc_pid_classique.json.
% Duree : deux a quatre minutes (auto-tests compris). Compatible R2024a.
% N'utilise pas la Fuzzy Logic Toolbox.
% ---------------------------------------------------------------------

clear; clc;

SOURCE = 'Buck_Commun';                           % modele commun (jamais modifie)
MDL    = 'Buck_Commun_Fuzzy_PID';                 % modele construit
TC_TXT = '1/(22000*10)';                          % periode du regulateur (texte exact)
TC     = 1/(22000*10);                            % la meme, en nombre (s)
CLASSE = 'ordonnanceur_flou_zhao';                % classe du bloc MATLAB System
NOM_SYS   = 'Ordonnanceur Flou Zhao';
NOM_ZOHE  = 'Echantillonneur Erreur Fuzzy';
NOM_DEMUX = 'Demux Gains Fuzzy';
NOM_N     = 'Constante N Filtre';
PWM       = sprintf('PWM Generator\n(DC-DC)');    % vrai saut de ligne

% --- Prerequis ---
for f = {[SOURCE '.slx'], 'charger_scenario.m', 'scenario_S1.mat', 'scenario_S8b.mat', [CLASSE '.m'], ...
         'fuzzy_pid_reglages.mat', 'predictions_banc_fuzzy_pid.json', 'predictions_banc_pid_classique.json'}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
    if numel(which(f{1}, '-all')) > 1
        error('Il doit exister un seul %s sur le chemin MATLAB : %s', f{1}, strjoin(which(f{1}, '-all'), ', '));
    end
end
if bdIsLoaded(SOURCE) && strcmp(get_param(SOURCE, 'Dirty'), 'on')
    error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', SOURCE);
end
for m = {MDL, SOURCE}
    if bdIsLoaded(m{1})
        close_system(m{1}, 0);                    % fermer sans enregistrer
    end
end
clear(CLASSE);                                    % recharger la classe si elle a change
charger_scenario('S1');                           % variables sc_* presentes avant le chargement

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

%% 2. Le PID, sa source et sa destination
pid = [MDL '/PID Controller'];
for b = {pid, [MDL '/Sum1'], [MDL '/' PWM]}
    if getSimulinkBlockHandle(b{1}) == -1
        close_system(MDL, 0);
        error('Bloc introuvable : %s', strrep(b{1}, newline, ' '));
    end
end
if ~strcmp(get_param(pid, 'ControllerParametersSource'), 'internal') || ~strcmp(get_param(pid, 'Form'), 'Parallel')
    close_system(MDL, 0);
    error('Le PID de Buck_Commun.slx devrait avoir des parametres internes et la forme Parallel.');
end
lh = get_param(pid, 'LineHandles');
if numel(lh.Inport) ~= 1
    close_system(MDL, 0);
    error('Le PID de Buck_Commun.slx devrait avoir une seule entree.');
end
[src_err, ~] = extremites_ligne(lh.Inport(1));    % ce qui alimente le PID (Sum1, sortie 1)
[~, dst_cmd] = extremites_ligne(lh.Outport(1));   % ce que le PID alimente (entree du PWM)
if ~strcmp(src_err.bloc, [MDL '/Sum1']) || numel(dst_cmd) ~= 1 || ~strcmp(dst_cmd(1).bloc, [MDL '/' PWM])
    close_system(MDL, 0);
    error('Le PID devrait etre alimente par Sum1 et alimenter seulement le PWM.');
end
N_TXT = get_param(pid, 'N');                      % valeur du filtre, telle qu'ecrite dans le bloc
fprintf('PID de Ziegler-Nichols : P = %s, I = %s, D = %s, N = %s (ces champs ne serviront plus).\n', ...
        get_param(pid, 'P'), get_param(pid, 'I'), get_param(pid, 'D'), N_TXT);

%% 3. Parametres du PID en entrees externes
delete_line(lh.Inport(1));                        % seule ligne supprimee, refaite plus bas
set_param(pid, 'ControllerParametersSource', 'external');
ports = ports_par_nom(pid);                       % numero d'entree de u, P, I, D, N
for nom = {'u', 'P', 'I', 'D', 'N'}
    if ~isfield(ports, nom{1})
        noms = fieldnames(ports);
        close_system(MDL, 0);
        error('Entree "%s" introuvable sous le masque du PID ; entrees trouvees : %s.', nom{1}, strjoin(noms', ', '));
    end
end
ph_pid = get_param(pid, 'PortHandles');
if numel(ph_pid.Inport) ~= 5
    close_system(MDL, 0);
    error('Le PID a parametres externes devrait avoir 5 entrees (u, P, I, D, N), il en a %d.', numel(ph_pid.Inport));
end
fprintf('PID en parametres externes : entrees u = %d, P = %d, I = %d, D = %d, N = %d.\n', ...
        ports.u, ports.P, ports.I, ports.D, ports.N);

%% 4. Ordonnanceur flou et branchements
pos = get_param(pid, 'Position');                 % [gauche haut droite bas]
zohe  = [MDL '/' NOM_ZOHE];
sys   = [MDL '/' NOM_SYS];
demux = [MDL '/' NOM_DEMUX];
cn    = [MDL '/' NOM_N];
add_block('simulink/Discrete/Zero-Order Hold', zohe, 'Position', [pos(1)-260, pos(4)+40, pos(1)-220, pos(4)+70]);
set_param(zohe, 'SampleTime', TC_TXT);
add_block('simulink/User-Defined Functions/MATLAB System', sys, 'Position', [pos(1)-190, pos(4)+30, pos(1)-110, pos(4)+80]);
set_param(sys, 'System', CLASSE);
set_param(sys, 'SimulateUsing', 'Interpreted execution');
ph = get_param(sys, 'PortHandles');
if numel(ph.Inport) ~= 1 || numel(ph.Outport) ~= 1
    close_system(MDL, 0);
    error('Le bloc %s devrait avoir 1 entree et 1 sortie (il en a %d et %d).', NOM_SYS, numel(ph.Inport), numel(ph.Outport));
end
add_block('simulink/Signal Routing/Demux', demux, 'Position', [pos(1)-80, pos(4)+25, pos(1)-75, pos(4)+85]);
set_param(demux, 'Outputs', '3');
add_block('simulink/Sources/Constant', cn, 'Position', [pos(1)-80, pos(4)+100, pos(1)-40, pos(4)+120]);
set_param(cn, 'Value', N_TXT);

relier(MDL, src_err.bloc, src_err.port, pid, ports.u);   % Sum1 -> entree u du PID (refaite)
relier(MDL, src_err.bloc, src_err.port, zohe, 1);        % Sum1 -> echantillonneur (branche)
relier(MDL, zohe, 1, sys, 1);
relier(MDL, sys, 1, demux, 1);
relier(MDL, demux, 1, pid, ports.P);
relier(MDL, demux, 2, pid, ports.I);
relier(MDL, demux, 3, pid, ports.D);
relier(MDL, cn, 1, pid, ports.N);
fprintf('Ordonnanceur "%s" (%s, execution interpretee) branche sur P, I, D ; N = %s constant.\n', NOM_SYS, CLASSE, N_TXT);

%% 5. Hygiene et sauvegarde
set_param(MDL, 'ReturnWorkspaceOutputs', 'off');
save_system(MDL);
close_system(MDL, 0);
fprintf('%s.slx sauvegarde.\n', MDL);

%% 6. Auto-tests
fprintf('\n=== Auto-tests ===\n');
pred_zn = jsondecode(fileread(fullfile(pwd, 'predictions_banc_pid_classique.json')));
pred = jsondecode(fileread(fullfile(pwd, 'predictions_banc_fuzzy_pid.json')));
reg = load(fullfile(pwd, 'fuzzy_pid_reglages.mat'));
pas_1ms = round(1e-3 / TC);                       % 220 periodes Tc
for essai = {{'T0', 'S1', 0.020, 0.1, true}, {'T1', 'S1', 0.020, 0.5, false}, {'T2', 'S8b', 0.020, 0.5, false}}
    nom = essai{1}{1}; code = essai{1}{2}; duree = essai{1}{3}; tol = essai{1}{4}; fixes = essai{1}{5};
    charger_scenario(code);
    r = simuler_court(MDL, NOM_SYS, duree, fixes);
    if fixes
        vb = pred_zn.essais.(code).v_toutes_les_ms(:);
        quoi = 'PID classique du banc';
    else
        vb = pred.essais.(code).v_toutes_les_ms(:);
        quoi = 'Fuzzy-PID du banc';
    end
    vs = r.v(1:pas_1ms:end);
    n = min(numel(vs), numel(vb));
    ecart = max(abs(vs(1:n) - vb(1:n)));
    Ks = r.K(1 + (n - 1) * pas_1ms, :);           % gains a (n-1) ms
    fprintf('  %s (%s, %.0f ms, %s) : ecart maximal sur Vout (toutes les ms) %.1f mV ; gains a %d ms : x[%.3f %.3f %.3f] de ZN', ...
            nom, code, duree * 1e3, quoi, ecart * 1e3, n - 1, Ks(1) / reg.K_ZN(1), Ks(2) / reg.K_ZN(2), Ks(3) / reg.K_ZN(3));
    if ~fixes
        Kb = pred.essais.(code).K_toutes_les_ms(n, :);
        fprintf(' (banc x[%.3f %.3f %.3f])', Kb(1) / reg.K_ZN(1), Kb(2) / reg.K_ZN(2), Kb(3) / reg.K_ZN(3));
    end
    fprintf('\n');
    if ~(ecart <= tol)
        if fixes
            error(['%s : ecart de %.3f V au PID classique avec les gains de Ziegler-Nichols : le branchement ' ...
                   'des entrees externes du PID est faux. Envoie cette sortie.'], nom, ecart);
        end
        error(['%s : ecart au banc de %.3f V, au-dela de %.2f V. Si T0 est passe, le bloc PID ne traite pas des ' ...
               'gains variables comme le banc le suppose : envoie cette sortie.'], nom, ecart, tol);
    end
    if strcmp(nom, 'T1')                          % retour au banc apres le demarrage
        ecart_fin = max(abs(vs(11:n) - vb(11:n)));    % de 10 ms a la fin (instants 10, 11, ... ms)
        fprintf('       de 10 a %d ms : ecart maximal %.1f mV (0.01 V au plus)\n', n - 1, ecart_fin * 1e3);
        if ~(ecart_fin <= 0.01)
            error('T1 : de 10 a %d ms, ecart au banc de %.3f V, au-dela de 0.01 V : envoie cette sortie.', n - 1, ecart_fin);
        end
    end
    if fixes && max(abs(r.K(end, :) - reg.K_ZN(:)') ./ reg.K_ZN(:)') > 0
        error('T0 : l''ordonnanceur n''a pas sorti les gains de Ziegler-Nichols.');
    end
end
fprintf(['\n%s.slx est construit et ses auto-tests sont passes.\n' ...
         'Lancer ensuite Simuler_Fuzzy_PID.m.\n'], MDL);


%% ===================== Fonctions =====================================

function ports = ports_par_nom(bloc)
    % Numero de chaque entree du bloc, d'apres le nom des Inport sous son
    % masque (lien de bibliotheque suivi).
    ib = find_system(bloc, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'SearchDepth', 1, 'BlockType', 'Inport');
    ports = struct();
    for i = 1:numel(ib)
        nom = get_param(ib{i}, 'Name');
        if isvarname(nom)
            ports.(nom) = str2double(get_param(ib{i}, 'Port'));
        end
    end
end

function r = simuler_court(mdl, nom_sys, duree, gains_fixes)
    % Simulation avec enregistrements temporaires (Vout et gains), modele
    % ferme SANS enregistrer, meme en cas d'erreur.
    load_system(fullfile(pwd, [mdl '.slx']));
    if gains_fixes                                % propriete logique : 'true' ou 'on' selon la release de MATLAB
        try
            set_param([mdl '/' nom_sys], 'GainsFixes', 'true');
        catch
            set_param([mdl '/' nom_sys], 'GainsFixes', 'on');
        end
    end
    ajouter_enregistrement(mdl, [mdl '/Log Vout tmp'], 'tmp_v', [mdl '/Vout'], 1);
    ajouter_enregistrement(mdl, [mdl '/Log K tmp'], 'tmp_K', [mdl '/' nom_sys], 1);
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
    r.v = lire(sortie, 'tmp_v');
    K = lire(sortie, 'tmp_K');
    if size(K, 2) ~= 3                            % composantes en lignes : on transpose
        K = K.';
    end
    r.K = K;
end

function M = lire(sortie, nom)
    % Lit un enregistrement To Workspace (scalaire ou vecteur) dans l'objet renvoye par sim().
    if ~any(strcmp(sortie.who, nom))
        error('Enregistrement "%s" absent des resultats de simulation.', nom);
    end
    donnee = sortie.get(nom);
    if isa(donnee, 'timeseries')
        M = double(squeeze(donnee.Data));
    else
        M = double(squeeze(donnee.signals.values));
    end
    if isvector(M)
        M = M(:);
    end
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
