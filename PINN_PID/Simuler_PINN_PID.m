%% SIMULER_PINN_PID.m
%
% Simulation du PINN-PID (modele Buck_Commun_PINN_PID.slx) sur les essais
% communs. Reglages propres a la methode : modele, bloc "PID Controller"
% dont la sortie est la commande, resultats attendus
% predictions_banc_pinn_pid.json et BLOC_GAINS : les gains [P I D] sont
% enregistres a la sortie de "PINN-PID Adaptatif" et compares au banc. Les
% resultats sont ecrits dans resultats_Buck_Commun_PINN_PID.mat.
% Ordre : banc_pinn_pid.py ; Tester_PINN_PID_Rejeu.m ;
% Construction_PINN_PID.m ; ce script.
% Duree : vingt a quarante minutes (l'adaptateur tourne en execution
% interpretee, 44 001 appels par essai de 0.2 s, et 5 iterations d'Adam
% sur un horizon de 110 pas a chaque fenetre optimisee).
% Resultats attendus : tous les evenements reviennent dans la bande ; les
% gains sont optimises sur 7 a 49 fenetres de 0.5 ms par essai, pendant le
% demarrage et apres chaque evenement. Chaque decision depend du seuil de
% 0.1 V sur l'erreur efficace de la fenetre : un ecart de quelques mV entre
% Simulink et le banc peut faire apparaitre ou disparaitre une
% optimisation, et les gains suivent alors un autre chemin. La validation du
% bloc est le test de rejeu (Tester_PINN_PID_Rejeu.m). Attendus ici : memes
% retours dans la bande, memes ordres de grandeur des gains et des IAE,
% egalite point par point jusqu'a la premiere decision qui bascule.
%
% Essais simules : ceux de scenarios_communs.json (onze essais de
% developpement), puis les essais complementaires (S10, hors du cout J), avec les memes grandeurs et la meme
% comparaison au banc. Reglage ESSAIS_A_SIMULER en tete du script : {} pour
% tous les essais, {'S10'} pour S10 seul (les resultats vont alors dans
% resultats_<MODELE>_S10.mat).
%
% OBJECTIF
% --------
% Simuler un modele construit sur la base commune (Buck_Commun.slx, ou le
% modele d'une methode construit a partir de lui) sur les onze essais S1 a
% S9, mesurer les grandeurs communes et, si un fichier de resultats
% attendus du banc Python est fourni, les comparer. Ce script sert a
% toutes les methodes : seuls les reglages du debut changent.
%
% COMMENT LES SIGNAUX SONT ENREGISTRES
% -------------------------------------
% Pour chaque essai : charger_scenario ecrit les profils dans le workspace
% de base ; des blocs To Workspace temporaires (Vout, sortie du regulateur,
% signal qui entre dans le PWM, courant iL, et la deuxieme sortie du
% regulateur s'il en a une) sont ajoutes en memoire ; le
% modele est simule puis ferme SANS etre enregistre. Les fichiers .slx
% verifies ne sont jamais modifies.
% Les blocs To Workspace propres au modele sont mis en commentaire pendant
% la simulation, en memoire seulement : ceux de Buck_Commun.slx
% heritent du pas du powergui, soit 26.4 millions de points par seconde
% simulee et par bloc, qui ralentiraient chaque essai sans servir ici.
%
% LES GRANDEURS (memes definitions que banc_commun.py, fonction grandeurs)
% -------------------------------------------------------------------------
%   e = consigne - Vout (vraie tension, pas la mesure bruitee).
%   Demarrage (0 a min(premier evenement, 50 ms)) : temps de montee de 10 a
%   90 V, depassement (% de 100 V), temps d'etablissement dans +-1 V,
%   erreur moyenne et ondulation crete a crete de 30 ms a la fin du
%   demarrage.
%   De 30 ms a la fin : IAE, ISE, ITAE, ecart maximal, ecart efficace.
%   Pour chaque evenement : IAE, ecart maximal, temps de retour dans +-1 V
%   et, sur les 10 dernieres ms avant l'evenement suivant, amplitude crete a
%   crete de l'ecart e et retour ou non dans la bande. Si l'ecart depasse
%   encore 1 V dans ces 10 dernieres ms, la tension n'est "pas revenue" et
%   le temps de retour vaut NaN (le banc ecrit null).
%   Activite de la commande (30 ms a la fin) : moyenne de |d(k) - d(k-1)|,
%   part du temps en butee, ecart type de d ; d est le signal qui entre
%   dans le PWM, qui doit rester dans [0.01 ; 0.99] (le script s'arrete
%   sinon).
%   Conduction discontinue : part du temps ou iL < 0.01 A.
%   Courant de la bobine iL (observe seulement, il n'entre dans aucun
%   critere) : courant crete pendant le demarrage, ondulation crete a crete
%   de 30 ms a la fin du demarrage, courant maximal et minimal dans la
%   fenetre de chaque evenement.
%   Temps de calcul : duree de la simulation Simulink divisee par le nombre
%   de periodes du regulateur (44 001 pour un essai de 0.2 s), et non par
%   le nombre de pas du powergui (120 fois plus).
%
% CE QUE PRODUIT CE SCRIPT
% -------------------------
%   - dans la console, un resume par essai et la comparaison au banc :
%     ecart sur Vout et sur iL echantillonnes toutes les ms, IAE, ecart
%     efficace, courant crete, gains (si le regulateur les sort et si le
%     banc les donne) et, pour chaque evenement, IAE, amplitude en fin de
%     fenetre et retour.
%   - resultats_<MODELE>.mat : signaux et grandeurs de chaque essai (et
%     gains du regulateur s'il les sort), pour la comparaison finale entre
%     methodes.
%
% Prerequis dans le dossier courant MATLAB : le modele, charger_scenario.m,
% scenarios_communs.json et les fichiers scenario_*.mat ; facultatif : le
% fichier de resultats attendus du banc.
% Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

clear; clc;

% --- Les trois reglages propres a chaque methode ---
MODELE = 'Buck_Commun_PINN_PID';                           % modele a simuler
BLOC_REGULATEUR = 'PID Controller';                       % bloc dont la sortie 1 est la commande
FICHIER_PREDICTIONS = 'predictions_banc_pinn_pid.json';    % '' : pas de comparaison
BLOC_GAINS = 'PINN-PID Adaptatif';                         % sortie 1 : gains [P I D] ('' : pas de gains a part)

TE = 1/(22000*10);                                        % periode commune des enregistrements (s)
TE_TXT = '1/(22000*10)';                                  % la meme, en texte exact pour Simulink
FICHIER_RESULTATS = ['resultats_' MODELE '.mat'];

% --- Essais a simuler ---
% {} : tous les essais de scenarios_communs.json, les onze essais du cout J
% puis les essais complementaires (S10, hors du cout J).
% Sinon, une liste de codes : {'S10'} ne simule que S10, sans relancer les
% onze essais deja valides. Avec une liste, les resultats sont ecrits dans
% resultats_<MODELE>_<codes>.mat (le fichier des onze essais n'est pas
% ecrase).
ESSAIS_A_SIMULER = {};

% --- Prerequis, un seul exemplaire de chaque fichier sur le chemin ---
for f = {[MODELE '.slx'], 'charger_scenario.m'}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
    copies = which(f{1}, '-all');
    if numel(copies) ~= 1
        error('Il doit exister exactement un %s sur le chemin MATLAB, il y en a %d :\n%s', ...
              f{1}, numel(copies), strjoin(copies, newline));
    end
end
if ~isfile(fullfile(pwd, 'scenarios_communs.json'))
    error('scenarios_communs.json introuvable (lancer scenarios_communs.py).');
end
liste = jsondecode(fileread(fullfile(pwd, 'scenarios_communs.json')));
if ~isfield(liste, 'format') || liste.format ~= 2
    error('scenarios_communs.json n''a pas le format attendu (voir scenarios_communs.py).');
end
if iscell(liste.essais)                                    % jsondecode rend une cellule ou une structure
    codes = cellfun(@(x) x.code, liste.essais, 'UniformOutput', false);
else
    codes = {liste.essais.code};
end
codes = reshape(codes, 1, []);                             % les onze essais du cout J
if isfield(liste, 'essais_complementaires') && ~isempty(liste.essais_complementaires)
    ec = liste.essais_complementaires;                     % essais hors du cout J (S10)
    if iscell(ec)
        codes = [codes, reshape(cellfun(@(x) x.code, ec, 'UniformOutput', false), 1, [])];
    else
        codes = [codes, {ec.code}];
    end
end
if ~isempty(ESSAIS_A_SIMULER)
    ESSAIS_A_SIMULER = cellstr(ESSAIS_A_SIMULER);
    inconnus = setdiff(ESSAIS_A_SIMULER, codes);
    if ~isempty(inconnus)
        error('ESSAIS_A_SIMULER : essai(s) inconnu(s) : %s (essais disponibles : %s).', ...
              strjoin(inconnus, ', '), strjoin(codes, ', '));
    end
    codes = codes(ismember(codes, ESSAIS_A_SIMULER));      % dans l'ordre de scenarios_communs.json
    FICHIER_RESULTATS = ['resultats_' MODELE '_' strjoin(codes, '_') '.mat'];
end
fprintf('Essais simules : %s (resultats dans %s).\n', strjoin(codes, ', '), FICHIER_RESULTATS);
pred = [];
if ~isempty(FICHIER_PREDICTIONS) && isfile(fullfile(pwd, FICHIER_PREDICTIONS))
    pred = jsondecode(fileread(fullfile(pwd, FICHIER_PREDICTIONS)));
    if ~isfield(pred, 'format') || pred.format ~= 2
        error(['%s n''a pas le format attendu : ' ...
               'relancer le banc Python de la methode.'], FICHIER_PREDICTIONS);
    end
    fprintf('Resultats attendus lus dans %s (%s).\n', FICHIER_PREDICTIONS, pred.regulateur);
else
    fprintf('Pas de fichier de resultats attendus : pas de comparaison au banc.\n');
end
if bdIsLoaded(MODELE)
    if strcmp(get_param(MODELE, 'Dirty'), 'on')
        error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', MODELE);
    end
    close_system(MODELE, 0);
end

% --- Simulation des essais ---
resultats = struct('code', {}, 'nom', {}, 't', {}, 'v', {}, 'consigne', {}, 'mu', {}, 'd', {}, ...
                   'iL', {}, 'K', {}, 'grandeurs', {}, 'duree_calcul_s', {}, 'us_par_pas', {});
for i = 1:numel(codes)
    sc = charger_scenario(codes{i});
    fprintf('\n=== %s : %s (%.0f ms) ===\n', sc.code, sc.nom, sc.duree * 1e3);
    load_system(fullfile(pwd, [MODELE '.slx']));
    bloc_reg = [MODELE '/' BLOC_REGULATEUR];
    if getSimulinkBlockHandle(bloc_reg) == -1
        close_system(MODELE, 0);
        error('Bloc regulateur introuvable : %s', bloc_reg);
    end
    src_pwm = source_entree_pwm(MODELE);
    ajouter_enregistrement(MODELE, [MODELE '/Log Vout sc'], 'sc_log_v', TE_TXT, [MODELE '/Vout'], 1);
    ajouter_enregistrement(MODELE, [MODELE '/Log mu sc'], 'sc_log_mu', TE_TXT, bloc_reg, 1);
    ajouter_enregistrement(MODELE, [MODELE '/Log d sc'], 'sc_log_d', TE_TXT, src_pwm.bloc, src_pwm.port);
    ajouter_enregistrement(MODELE, [MODELE '/Log iL sc'], 'sc_log_iL', TE_TXT, [MODELE '/iL'], 1);
    ph_reg = get_param(bloc_reg, 'PortHandles');
    if ~isempty(BLOC_GAINS)                                % gains sortis par un autre bloc (sortie 1)
        if getSimulinkBlockHandle([MODELE '/' BLOC_GAINS]) == -1
            close_system(MODELE, 0);
            error('Bloc des gains introuvable : %s', BLOC_GAINS);
        end
        avec_gains = true;
        ajouter_enregistrement(MODELE, [MODELE '/Log K sc'], 'sc_log_K', TE_TXT, [MODELE '/' BLOC_GAINS], 1);
    else
        avec_gains = numel(ph_reg.Outport) >= 2;           % deuxieme sortie : gains d'une methode adaptative
        if avec_gains
            ajouter_enregistrement(MODELE, [MODELE '/Log K sc'], 'sc_log_K', TE_TXT, bloc_reg, 2);
        end
    end
    neutraliser_enregistrements_modele(MODELE);
    t_debut = tic;
    try
        % StopTime : un demi-pas apres le dernier instant du profil, pour que
        % cet instant soit enregistre quels que soient les arrondis ; les
        % enregistrements sont ensuite ramenes aux n instants du profil.
        sortie = sim(MODELE, 'StopTime', sprintf('%.17g', (numel(sc.t) - 0.5) * TE), 'ReturnWorkspaceOutputs', 'on');
    catch ME
        close_system(MODELE, 0);
        rethrow(ME);
    end
    duree_calcul = toc(t_debut);
    close_system(MODELE, 0);                              % SANS enregistrer

    [t, v]    = extraire_enregistrement(sortie, 'sc_log_v');
    [t2, mu]  = extraire_enregistrement(sortie, 'sc_log_mu');
    [t3, dpw] = extraire_enregistrement(sortie, 'sc_log_d');
    [t4, iL]  = extraire_enregistrement(sortie, 'sc_log_iL');
    K = [];
    if avec_gains
        [tk, K] = extraire_matrice(sortie, 'sc_log_K');   % une ligne par instant, une colonne par gain
        if size(K, 1) >= numel(sc.t)
            K = K(1:numel(sc.t), :);
            tk = tk(1:numel(sc.t));
        end
        if size(K, 1) ~= numel(sc.t) || max(abs(tk - sc.t)) > 1e-9
            error('%s : enregistrement des gains non aligne avec le profil.', sc.code);
        end
    end
    n_sc = numel(sc.t);                                   % instants du profil
    if numel(t) >= n_sc && numel(t2) >= n_sc && numel(t3) >= n_sc && numel(t4) >= n_sc
        t = t(1:n_sc); t2 = t2(1:n_sc); t3 = t3(1:n_sc); t4 = t4(1:n_sc);
        v = v(1:n_sc); mu = mu(1:n_sc); dpw = dpw(1:n_sc); iL = iL(1:n_sc);
    end
    if ~isequal(numel(t), numel(t2), numel(t3), numel(t4), numel(sc.t)) || ...
            max(abs([t - t2; t - t3; t - t4; t - sc.t])) > 1e-9
        error('%s : enregistrements non alignes avec le profil (%d pas attendus, %d obtenus).', ...
              sc.code, numel(sc.t), numel(t));
    end
    consigne = 100 + sc.dvref;
    % Regle commune : une seule saturation physique [0.01 ; 0.99], juste avant
    % le PWM. Si le signal qui entre dans le PWM en sort, le modele ne suit
    % pas le protocole (et le banc, qui borne a [0.01 ; 0.99], ne le
    % reproduirait pas).
    if any(dpw < 0.01 - 1e-9 | dpw > 0.99 + 1e-9)
        error(['%s : le signal qui entre dans le PWM sort de [0.01 ; 0.99] (de %.4f a %.4f). ' ...
               'Il manque la saturation physique devant le PWM.'], sc.code, min(dpw), max(dpw));
    end
    d = dpw;
    g = grandeurs(v, consigne, d, iL, sc.evenements, TE);
    resultats(end + 1) = struct('code', sc.code, 'nom', sc.nom, 't', t, 'v', v, 'consigne', consigne, ...
                                'mu', mu, 'd', d, 'iL', iL, 'K', K, 'grandeurs', g, 'duree_calcul_s', duree_calcul, ...
                                'us_par_pas', duree_calcul / numel(t) * 1e6); %#ok<SAGROW>
    fprintf('  Simulink : %s\n', resume(g));
    fprintf('  temps de calcul : %.1f s (%.1f us par periode du regulateur)\n', duree_calcul, duree_calcul / numel(t) * 1e6);
    for j = 1:numel(g.evenements)
        ev = g.evenements(j);
        if ev.reste_dans_bande
            etat = 'reste dans la bande';
        elseif ev.revenu
            etat = sprintf('retour en %.2f ms', ev.t_retour_ms);
        else
            etat = 'PAS REVENU dans la bande avant l''evenement suivant';
        end
        fprintf('    evenement a %6.1f ms : IAE %7.2f mV.s, ecart max %6.2f V, ecart crete a crete en fin de fenetre %5.2f V, %s ; iL de %.1f a %.1f A\n', ...
                ev.t_ms, ev.IAE * 1e3, ev.e_max_V, ev.e_crete_a_crete_fin_V, etat, ev.iL_min_A, ev.iL_max_A);
    end
    if ~isempty(pred) && isfield(pred.essais, sc.code)
        gb = pred.essais.(sc.code).grandeurs;
        vb = pred.essais.(sc.code).v_toutes_les_ms(:);
        vs = v(1:round(1e-3 / TE):end);
        n_cmp = min(numel(vb), numel(vs));
        fprintf('  banc     : %s\n', resume(gb));
        fprintf('  ecart Simulink - banc sur Vout (toutes les ms) : %.1f mV efficace, %.1f mV au plus\n', ...
                sqrt(mean((vs(1:n_cmp) - vb(1:n_cmp)).^2)) * 1e3, max(abs(vs(1:n_cmp) - vb(1:n_cmp))) * 1e3);
        fprintf('  IAE : Simulink %.2f, banc %.2f mV.s (%+.1f %%) ; ecart efficace : Simulink %.1f, banc %.1f mV (%+.1f %%)\n', ...
                g.IAE * 1e3, gb.IAE * 1e3, (g.IAE / gb.IAE - 1) * 100, ...
                g.e_eff_V * 1e3, gb.e_eff_V * 1e3, (g.e_eff_V / gb.e_eff_V - 1) * 100);
        if isfield(pred.essais.(sc.code), 'iL_toutes_les_ms') && isfield(gb, 'iL_max_dem_A')
            ib = pred.essais.(sc.code).iL_toutes_les_ms(:);     % courant du banc toutes les ms
            is = iL(1:round(1e-3 / TE):end);
            n_i = min(numel(ib), numel(is));
            fprintf(['  courant iL (toutes les ms) : ecart Simulink - banc %.1f mA efficace, %.1f mA au plus ; ' ...
                     'crete au demarrage : Simulink %.2f A, banc %.2f A\n'], ...
                    sqrt(mean((is(1:n_i) - ib(1:n_i)).^2)) * 1e3, max(abs(is(1:n_i) - ib(1:n_i))) * 1e3, ...
                    g.iL_max_dem_A, gb.iL_max_dem_A);
        end
        if ~isempty(K) && isfield(pred.essais.(sc.code), 'K_toutes_les_ms')
            Kb = pred.essais.(sc.code).K_toutes_les_ms;          % gains du banc toutes les ms
            Ks = K(1:round(1e-3 / TE):end, :);
            n_k = min(size(Kb, 1), size(Ks, 1));
            ecart_K = max(max(abs(Ks(1:n_k, :) - Kb(1:n_k, :)) ./ abs(Kb(1:n_k, :))));
            fprintf(['  gains : ecart relatif maximal Simulink - banc (toutes les ms) %.2g ; gains finaux ' ...
                     'Simulink [%.5f %.4e %.4e], banc [%.5f %.4e %.4e]\n'], ecart_K, K(end, 1), K(end, 2), ...
                    K(end, 3), Kb(end, 1), Kb(end, 2), Kb(end, 3));
        end
        evb = gb.evenements;                              % structure, cellule ou vide selon jsondecode
        if iscell(evb)
            evb = [evb{:}];
        end
        oscille = false;                                  % le banc predit-il une oscillation entretenue ?
        for j = 1:min(numel(evb), numel(g.evenements))
            fprintf('    evenement a %6.1f ms : IAE %7.2f / %7.2f mV.s ; ecart crete a crete en fin %5.2f / %5.2f V ; retour : %s / %s (Simulink / banc)\n', ...
                    g.evenements(j).t_ms, g.evenements(j).IAE * 1e3, evb(j).IAE * 1e3, ...
                    g.evenements(j).e_crete_a_crete_fin_V, evb(j).e_crete_a_crete_fin_V, ...
                    texte_retour(g.evenements(j)), texte_retour(evb(j)));
            oscille = oscille || ~evb(j).revenu || evb(j).e_crete_a_crete_fin_V > 0.5;
        end
        if oscille
            fprintf(['  NB : le banc predit ici une oscillation entretenue du regulateur. L''ecart point par point\n' ...
                     '  depend alors de la phase de l''oscillation : comparer plutot IAE, amplitudes et retours.\n']);
        end
        if sc.q > 1e-6
            fprintf(['  NB : mesure quantifiee. Un ecart d''arrondi de l''ordre de 1e-15 V suffit a changer une\n' ...
                     '  decision du CAN : les trajectoires se separent point par point, seuls IAE et ecart\n' ...
                     '  efficace se comparent.\n']);
        end
    end
end

save(fullfile(pwd, FICHIER_RESULTATS), 'resultats', 'MODELE', 'BLOC_REGULATEUR', '-v7.3');
fprintf('\n%s ecrit.\n', FICHIER_RESULTATS);

% --- Tableau recapitulatif ---
fprintf('\n==========================================================================================\n');
fprintf(' %s : recapitulatif (Simulink)\n', MODELE);
fprintf('==========================================================================================\n');
fprintf('  %-5s %-24s %9s %9s %10s %9s %9s %9s %7s %9s\n', 'essai', 'nom', 'depass.%', 'etabli ms', ...
        'IAE mV.s', 'e max V', 'e eff mV', '|dd|', 'DCM %', 'iL dem A');
for i = 1:numel(resultats)
    g = resultats(i).grandeurs;
    fprintf('  %-5s %-24s %9.2f %9.2f %10.2f %9.2f %9.1f %9.4f %7.1f %9.2f\n', resultats(i).code, resultats(i).nom, ...
            g.depassement_pct, g.t_etab_ms, g.IAE * 1e3, g.e_max_V, g.e_eff_V * 1e3, g.dd_moyen, g.dcm_pct, ...
            g.iL_max_dem_A);
end


%% ===================== Fonctions =====================================

function o = grandeurs(v, consigne, d, iL, evenements, Te)
    % Memes definitions que la fonction grandeurs de banc_commun.py. Les
    % numeros de pas commencent a 0 (comme en Python) ; l'indice MATLAB du
    % pas n est n + 1.
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
    if isfield(g, 'iL_max_dem_A')                         % fichiers de predictions anterieurs : sans courant
        s = sprintf('%s ; iL crete au demarrage %.1f A', s, g.iL_max_dem_A);
    end
end

function s = texte_retour(ev)
    % Etat de retour dans la bande d'un evenement (Simulink ou banc).
    if ev.reste_dans_bande
        s = 'dans la bande';
    elseif ev.revenu
        s = sprintf('%.2f ms', ev.t_retour_ms);
    else
        s = 'pas revenu';
    end
end

function s = source_entree_pwm(mdl)
    % Bloc et numero de port qui alimentent l'entree du PWM Generator.
    pwm = [mdl '/' sprintf('PWM Generator\n(DC-DC)')];  % vrai saut de ligne
    if getSimulinkBlockHandle(pwm) == -1
        error('PWM Generator (DC-DC) introuvable dans %s.', mdl);
    end
    ph = get_param(pwm, 'PortHandles');
    ligne = get_param(ph.Inport(1), 'Line');
    if ligne == -1
        error('L''entree du PWM n''est reliee a rien dans %s.', mdl);
    end
    s.bloc = getfullname(get_param(ligne, 'SrcBlockHandle'));
    s.port = get_param(get_param(ligne, 'SrcPortHandle'), 'PortNumber');
end

function neutraliser_enregistrements_modele(mdl)
    % Met en commentaire, en memoire seulement, les blocs To Workspace de
    % premier niveau qui ne sont pas les enregistrements de ce script
    % (variables 'sc_log_...'). Le modele est ferme sans enregistrer.
    blocs = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
    for i = 1:numel(blocs)
        if ~startsWith(get_param(blocs{i}, 'VariableName'), 'sc_log_')
            set_param(blocs{i}, 'Commented', 'on');
        end
    end
end

function relier(mdl, bloc_src, port_src, bloc_dst, port_dst)
    % Ligne de signal entre deux ports (handles demandes au moment meme).
    ph_src = get_param(bloc_src, 'PortHandles');
    ph_dst = get_param(bloc_dst, 'PortHandles');
    add_line(mdl, ph_src.Outport(port_src), ph_dst.Inport(port_dst), 'autorouting', 'on');
end

function ajouter_enregistrement(mdl, chemin, variable, Te_txt, bloc_src, port_src)
    % Bloc To Workspace temporaire (MaxDataPoints = inf), echantillonne a
    % Te_txt (texte exact).
    if getSimulinkBlockHandle(chemin) ~= -1
        delete_block(chemin);
    end
    pos = get_param(bloc_src, 'Position');
    add_block('simulink/Sinks/To Workspace', chemin, 'Position', [pos(3)+120, pos(2)+40, pos(3)+200, pos(2)+65]);
    set_param(chemin, 'VariableName', variable, 'SaveFormat', 'Timeseries', ...
              'SampleTime', Te_txt, 'MaxDataPoints', 'inf');
    relier(mdl, bloc_src, port_src, chemin, 1);
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
