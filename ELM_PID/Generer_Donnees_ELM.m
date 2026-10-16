%% GENERER_DONNEES_ELM.m
%
% OBJECTIF
% --------
% Enregistrer le vrai circuit Simscape en boucle ouverte pour
% l'apprentissage et la validation de l'ELM : simuler Donnees_ELM.slx
% (construit par Construction_Modele_Donnees_ELM.m) sur les quatre enregistrements prepares
% par excitation_donnees_elm.py, et comparer chacun au banc commun.
%
% LES QUATRE ENREGISTREMENTS (voir excitation_donnees_elm.py)
% ------------------------------------------------------------
%   ELM_A1, ELM_A2 : apprentissage, 1 s chacun (graines 31, 32) ;
%   ELM_V1, ELM_V2 : validation, 0.5 s chacun (graines 41, 42).
% Rapport cyclique en paliers de 0.5 a 20 ms, Vin de 150 a 240 V et charge
% de 4 a 98 ohms en rampes, tous les changements sur un debut de fenetre de
% 0.5 ms. L'ELM ne verra jamais les enregistrements de validation pendant
% son apprentissage.
%
% CE QUE PRODUIT CE SCRIPT
% -------------------------
% Un fichier donnees_elm_<code>.mat par enregistrement, contenant :
%   vout, il, d, vin : tension de sortie (V), courant de la bobine (A),
%     rapport cyclique applique au PWM et mesure de Vin (V), un echantillon
%     a chaque instant k*Tc (Tc = 1/220000 s), en simple precision
%     (7 chiffres significatifs, 0.01 mV a 100 V). La mesure vin a un pas
%     du powergui de retard sur le profil sc_vin ("Retard Vin") : au
%     premier instant d'une fenetre ou Vin change, elle lit encore la valeur
%     precedente. Pour Vin par fenetre, prendre le profil du fichier
%     scenario_<code>.mat ;
%   r_eq : resistance de charge equivalente a chaque instant (ohms), en
%     simple precision ;
%   Te, N, NF : periode, nombre d'instants, instants par fenetre (110) ;
%   code, role, graine, duree ; ecart_banc_fen_mV, ecart_banc_ms_mV ;
%   duree_calcul_s : duree du calcul.
% Le modele n'est jamais modifie sur le disque : les enregistrements sont
% ajoutes en memoire et le modele est ferme sans enregistrer.
%
% COMPARAISON AU BANC
% --------------------
% Pour chaque enregistrement : ecart entre Simulink et le banc commun sur
% la tension moyenne de chaque fenetre de 0.5 ms (ce que l'ELM apprend) et
% sur la tension prise toutes les millisecondes. Sur les essais communs, le
% banc a retrouve Simulink a moins de 1 mV efficace : on attend ici
% quelques mV au plus. Le script affiche ou se trouve l'ecart maximal, et
% l'ecart par classe de charge et dans les fenetres ou le courant de la
% bobine s'annule. Les donnees Simulink sont ecrites dans tous les cas.
%
% ORDRE D'EXECUTION (donnees de l'ELM-PID)
% -----------------------------------------
%   1. excitation_donnees_elm.py ; 2. Construction_Modele_Donnees_ELM.m ;
%   3. ce script.
%
% Prerequis dans le dossier courant MATLAB : Donnees_ELM.slx,
% charger_scenario.m, scenario_ELM_*.mat, predictions_banc_donnees_elm.mat.
% Duree : 3 s simulees au pas de 37.88 ns, soit 8 a 30 minutes selon la
% machine. Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

clear; clc;

MDL  = 'Donnees_ELM';
TC   = 1/(22000*10);                              % periode des enregistrements (s)
NF   = 110;                                       % instants par fenetre de 0.5 ms
CODES = {'ELM_A1', 'ELM_A2', 'ELM_V1', 'ELM_V2'};
ROLES = {'apprentissage', 'apprentissage', 'validation', 'validation'};

% --- Prerequis ---
chemin_mdl = fullfile(pwd, [MDL '.slx']);
if ~isfile(chemin_mdl)
    error('%s.slx introuvable : lance d''abord Construction_Modele_Donnees_ELM.m.', MDL);
end
if numel(which([MDL '.slx'], '-all')) ~= 1
    error('Il doit exister exactement un %s.slx sur le chemin MATLAB.', MDL);
end
for i = 1:numel(CODES)
    if ~isfile(fullfile(pwd, ['scenario_' CODES{i} '.mat']))
        error('scenario_%s.mat introuvable (excitation_donnees_elm.py).', CODES{i});
    end
end
P = load(fullfile(pwd, 'predictions_banc_donnees_elm.mat'));   % resultats attendus du banc
if bdIsLoaded(MDL)
    if strcmp(get_param(MDL, 'Dirty'), 'on')
        error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', MDL);
    end
    close_system(MDL, 0);
end

bilan = cell(numel(CODES), 1);
for i = 1:numel(CODES)
    code = CODES{i};
    sc = charger_scenario(code);                  % Vin, charge, R0 dans le workspace de base
    m = load(fullfile(pwd, ['scenario_' code '.mat']), 'dexc', 'graine');
    dexc = double(m.dexc(:));
    n = numel(sc.t);
    td = [0; ((1:n-1)' - 0.5) * TC];              % instants decales d'un demi-pas (charger_scenario)
    assignin('base', 'sc_dexc', [td, dexc]);
    fprintf('\n=== %s (%s, %.2f s, graine %d) ===\n', code, ROLES{i}, sc.duree, round(double(m.graine)));

    % Modele charge, enregistrements temporaires a Tc, journaux d'origine neutralises
    load_system(chemin_mdl);
    ajouter_enregistrement(MDL, [MDL '/Log Vout elm'], 'elm_log_v', [MDL '/Vout'], 1);
    ajouter_enregistrement(MDL, [MDL '/Log iL elm'], 'elm_log_iL', [MDL '/iL'], 1);
    ajouter_enregistrement(MDL, [MDL '/Log d elm'], 'elm_log_d', [MDL '/Commande Excitation'], 1);
    ajouter_enregistrement(MDL, [MDL '/Log Vin elm'], 'elm_log_vin', [MDL '/Vin'], 1);
    blocs = find_system(MDL, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
    for j = 1:numel(blocs)
        if ~startsWith(get_param(blocs{j}, 'VariableName'), 'elm_log_')
            set_param(blocs{j}, 'Commented', 'on');
        end
    end
    t0 = tic;
    try
        % StopTime : un demi-pas apres le dernier instant, pour qu'il soit enregistre
        % SaveTime 'off' : le vecteur des instants du solveur (26.4 millions de
        % points par seconde simulee) n'est pas garde ; les enregistrements
        % To Workspace ont leurs propres instants.
        sortie = sim(MDL, 'StopTime', sprintf('%.17g', (n - 0.5) * TC), 'ReturnWorkspaceOutputs', 'on', ...
                     'SaveTime', 'off');
    catch ME
        close_system(MDL, 0);
        rethrow(ME);
    end
    duree_calcul = toc(t0);
    close_system(MDL, 0);                         % SANS enregistrer

    [t, vout] = extraire_enregistrement(sortie, 'elm_log_v');
    [t2, il]  = extraire_enregistrement(sortie, 'elm_log_iL');
    [t3, d]   = extraire_enregistrement(sortie, 'elm_log_d');
    [t4, vin] = extraire_enregistrement(sortie, 'elm_log_vin');
    if numel(t) < n || numel(t2) < n || numel(t3) < n || numel(t4) < n
        error('%s : enregistrements trop courts (%d instants pour %d attendus).', code, numel(t), n);
    end
    t = t(1:n); vout = vout(1:n); il = il(1:n); d = d(1:n); vin = vin(1:n);
    if max(abs(t - sc.t)) > 1e-9 || max(abs([t2(1:n); t3(1:n); t4(1:n)] - [t; t; t])) > 1e-9
        error('%s : enregistrements non alignes avec les instants k*Tc.', code);
    end
    if max(abs(d - dexc)) > 1e-12
        error('%s : le rapport cyclique applique differe du profil impose (%.3g).', code, max(abs(d - dexc)));
    end

    % Comparaison au banc commun
    n_fen = floor(n / NF);
    v_fen = mean(reshape(vout(1:n_fen * NF), NF, n_fen), 1)';
    v_fen_banc = P.(['v_fen_' code])(:);
    v_ms_banc = P.(['v_ms_' code])(:);
    v_ms = vout(1:220:end);
    k = min(numel(v_ms), numel(v_ms_banc));
    e_fen = v_fen - v_fen_banc(1:n_fen);
    e_ms = v_ms(1:k) - v_ms_banc(1:k);
    fprintf('  temps de calcul %.0f s ; Vout de %.1f a %.1f V ; d de %.3f a %.3f\n', duree_calcul, ...
            min(vout), max(vout), min(d), max(d));
    fprintf('  ecart au banc : fenetres de 0.5 ms %.2f mV efficace, %.2f mV au plus ; toutes les ms %.2f mV efficace, %.2f mV au plus\n', ...
            sqrt(mean(e_fen.^2)) * 1e3, max(abs(e_fen)) * 1e3, sqrt(mean(e_ms.^2)) * 1e3, max(abs(e_ms)) * 1e3);

    % Diagnostic de l'ecart au banc. Les donnees Simulink sont ecrites dans
    % tous les cas : ce sont elles qui servent a l'apprentissage. L'ecart dit
    % seulement si le banc, qui servira a regler l'ELM-PID, reproduit le
    % circuit en boucle ouverte sur toute la plage.
    r_eq = 1 ./ (1 / sc.R0 + sc.gx);              % charge equivalente a chaque instant
    fen = @(x) mean(reshape(x(1:n_fen * NF), NF, n_fen), 1)';      % moyenne par fenetre
    r_fen = fen(r_eq);
    vin_fen = fen(sc.vin);                        % profil de Vin (sans le retard de la mesure)
    d_fen = fen(d);
    il_min = min(reshape(il(1:n_fen * NF), NF, n_fen), [], 1)';     % courant minimal de chaque fenetre
    [e_max, i_max] = max(abs(e_fen));
    fprintf('  ecart maximal a %.1f ms : %+.3f V (R = %.1f ohms, Vin = %.0f V, d = %.3f, Vout = %.1f V, iL min = %.3f A)\n', ...
            (i_max - 1) * NF * TC * 1e3, e_fen(i_max), r_fen(i_max), vin_fen(i_max), d_fen(i_max), v_fen(i_max), il_min(i_max));
    i_prem = find(abs(e_fen) > 0.05, 1);
    if ~isempty(i_prem)
        fprintf('  premiere fenetre au-dela de 0.05 V : %.1f ms (R = %.1f ohms, Vin = %.0f V, d = %.3f, iL min = %.3f A)\n', ...
                (i_prem - 1) * NF * TC * 1e3, r_fen(i_prem), vin_fen(i_prem), d_fen(i_prem), il_min(i_prem));
    end
    classes = {r_fen < 6, r_fen >= 6 & r_fen < 25, r_fen >= 25};
    noms_classes = {'R < 6 ohms', '6 a 25 ohms', 'R >= 25 ohms'};
    for c = 1:3
        if any(classes{c})
            fprintf('  %-13s : %4d fenetres, ecart efficace %7.1f mV\n', noms_classes{c}, sum(classes{c}), ...
                    sqrt(mean(e_fen(classes{c}).^2)) * 1e3);
        end
    end
    dcm = il_min < 0.01;                          % fenetres ou le courant s'annule
    e_dcm = 0;
    if any(dcm)
        e_dcm = sqrt(mean(e_fen(dcm).^2)) * 1e3;
    end
    fprintf('  fenetres avec iL < 10 mA : %d (ecart efficace %.1f mV) ; autres : ecart efficace %.1f mV\n', ...
            sum(dcm), e_dcm, sqrt(mean(e_fen(~dcm).^2)) * 1e3);
    if e_max > 0.05
        fprintf('  ATTENTION : ecart au banc au-dela de 0.05 V sur une fenetre ; fichier ecrit quand meme.\n');
    end

    % Ecriture (format v7, lisible par scipy.io.loadmat)
    donnees = struct('vout', single(vout), 'il', single(il), 'd', single(d), 'vin', single(vin), ...
                     'r_eq', single(r_eq), 'Te', TC, 'N', n, 'NF', NF, 'code', code, 'role', ROLES{i}, ...
                     'graine', double(m.graine), 'duree', sc.duree, ...
                     'ecart_banc_fen_mV', [sqrt(mean(e_fen.^2)), max(abs(e_fen))] * 1e3, ...
                     'ecart_banc_ms_mV', [sqrt(mean(e_ms.^2)), max(abs(e_ms))] * 1e3, ...
                     'ecart_banc_fen_V', single(e_fen), ...
                     'duree_calcul_s', duree_calcul);
    save(fullfile(pwd, ['donnees_elm_' code '.mat']), '-struct', 'donnees', '-v7');
    fprintf('  donnees_elm_%s.mat ecrit.\n', code);
    bilan{i} = sprintf('  %-7s %-14s %6.0f s   %7.2f / %7.2f mV   %7.2f / %7.2f mV', code, ROLES{i}, duree_calcul, ...
                       donnees.ecart_banc_fen_mV, donnees.ecart_banc_ms_mV);
end

fprintf('\n==========================================================================\n');
fprintf(' Donnees ELM : ecart au banc (efficace / maximal)\n');
fprintf('==========================================================================\n');
fprintf('  code    role            calcul   fenetres 0.5 ms         toutes les ms\n');
fprintf('%s\n', bilan{:});
fprintf(['\nAttendu : quelques mV au plus. Un ecart de l''ordre du volt signale un probleme de ' ...
         'construction ou de profil.\n']);


%% ===================== Fonctions =====================================

function relier(mdl, bloc_src, port_src, bloc_dst, port_dst)
    % Ligne de signal entre deux ports, handles demandes au moment meme.
    ph_src = get_param(bloc_src, 'PortHandles');
    ph_dst = get_param(bloc_dst, 'PortHandles');
    add_line(mdl, ph_src.Outport(port_src), ph_dst.Inport(port_dst), 'autorouting', 'on');
end

function ajouter_enregistrement(mdl, chemin, variable, bloc_src, port_src)
    % Bloc To Workspace temporaire a Tc (texte exact), MaxDataPoints = inf.
    if getSimulinkBlockHandle(chemin) ~= -1
        delete_block(chemin);
    end
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
