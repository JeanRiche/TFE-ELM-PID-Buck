%% TESTER_ELM_PID_REJEU.m
%
% Compare gains et decisions de fenetre du bloc MATLAB a ceux du banc : le
% bloc est l'adaptateur des gains du PID (entrees e, mesure, u ; sortie K).
% La commande u du banc est donnee au bloc.
%
% OBJECTIF
% --------
% Verifier, avant toute simulation Simulink, que le bloc MATLAB
% elm_pid_adaptatif.m calcule exactement ce que calcule la classe
% AdaptateurELM de banc_elm_pid.py. Le banc a enregistre, pas par pas,
% l'erreur e, la mesure et la commande u (sortie du bloc PID) vues par
% l'adaptateur sur cinq essais (S1, S3, S7b, S8b, S9), ainsi que les gains
% K qu'il a sortis. Ce script donne les memes e, mesure et u au bloc
% MATLAB, en dehors de Simulink, et compare les gains. Le bloc PID
% lui-meme est le bloc de la base commune : il est controle par l'auto-test
% T0 de Construction_ELM_PID.m (sans adaptation, PID classique).
%
% CE QUI EST COMPARE
% ------------------
%   - les vecteurs de test du modele (elm_pid_modele.mat) : prediction et
%     jacobien de evaluer_modele, contre les valeurs ecrites par
%     entrainement_elm.py ;
%   - pour chaque essai et chaque pas : les trois gains ;
%   - pour chaque fenetre de 0.5 ms : porte ouverte ou non, gains changes
%     ou non.
% Tolerances : 1e-9 en relatif sur les gains, aucune decision
% de fenetre differente. Les ecarts attendus viennent seulement de l'ordre
% des operations dans les produits matriciels (quelques 1e-15).
%
% ORDRE D'EXECUTION
% ------------------
%   1. entrainement_elm.py ; 2. ensemble_gains_elm.py ; 3. banc_elm_pid.py ;
%   4. ce script ; 5. Construction_ELM_PID.m ; 6. Simuler_ELM_PID.m.
%
% Prerequis dans le dossier courant MATLAB : elm_pid_adaptatif.m,
% elm_pid_modele.mat, elm_pid_reglages.mat, ensemble_gains_elm.mat,
% reference_rejeu_elm_pid.mat.
% Duree : une demi-minute a deux minutes (230 000 pas). Compatible R2024a.
% ---------------------------------------------------------------------

clear; clc;

% --- Prerequis ---
for f = {'elm_pid_adaptatif.m', 'elm_pid_modele.mat', 'elm_pid_reglages.mat', 'ensemble_gains_elm.mat', ...
         'reference_rejeu_elm_pid.mat'}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
    if numel(which(f{1}, '-all')) ~= 1
        error('Il doit exister exactement un %s sur le chemin MATLAB.', f{1});
    end
end
clear elm_pid_adaptatif                           % recharger la classe si elle a change

% --- 1. Vecteurs de test du modele ---
m = load('elm_pid_modele.mat');
ecart_y = 0; ecart_J = 0;
for i = 1:size(m.X_TEST, 1)
    [y, J] = elm_pid_adaptatif.evaluer_modele(double(m.W), double(m.B(:)'), double(m.BETA(:)), ...
        double(m.X_MOY(:)'), double(m.X_EC(:)'), double(m.T_MOY), double(m.T_EC), double(m.X_TEST(i, :)));
    ecart_y = max(ecart_y, abs(y - m.Y_TEST(i)));
    ecart_J = max(ecart_J, abs(J - m.J_TEST(i, 1)));
end
fprintf('Vecteurs de test du modele : ecart %.1e V sur la prediction, %.1e sur le jacobien.\n', ecart_y, ecart_J);
if ecart_y > 1e-9 || ecart_J > 1e-9
    error('evaluer_modele ne redonne pas les vecteurs de test de entrainement_elm.py.');
end

% --- 2. Rejeu pas par pas ---
ref = load('reference_rejeu_elm_pid.mat');
codes = {'S1', 'S3', 'S7b', 'S8b', 'S9'};
tout_ok = true;
for c = 1:numel(codes)
    r = ref.(codes{c});
    e = double(r.e(:));
    mes = double(r.mesure(:));
    ub = double(r.u(:));
    N = numel(e);
    bloc = elm_pid_adaptatif('AfficherBilan', false);
    K = zeros(N, 3);
    porte = false(0, 1);
    adapte = false(0, 1);
    nfen = 0;
    t0 = tic;
    for k = 1:N
        K(k, :) = bloc(e(k), mes(k), ub(k));     % gains du pas, puis mise a jour
        info = lire_fenetre(bloc);
        if info.nfen > nfen                       % une fenetre vient de se terminer
            nfen = info.nfen;
            porte(end + 1, 1) = info.porte; %#ok<SAGROW>
            adapte(end + 1, 1) = info.adapte; %#ok<SAGROW>
        end
    end
    duree = toc(t0);
    release(bloc);
    ecart_K = max(max(abs(K - double(r.K)) ./ abs(double(r.K))));
    nw = min(numel(porte), numel(r.porte));
    diff_porte = sum(porte(1:nw) ~= logical(r.porte(1:nw)));
    diff_adapte = sum(adapte(1:nw) ~= logical(r.adapte(1:nw)));
    ok = ecart_K <= 1e-9 && diff_porte == 0 && diff_adapte == 0 && numel(porte) == numel(r.porte);
    tout_ok = tout_ok && ok;
    fprintf(['%-4s : %6d pas, %4d fenetres ; ecart max sur les gains %.1e (relatif) ; ' ...
             'portes differentes %d, adaptations differentes %d ; %4.1f s  -> %s\n'], codes{c}, N, numel(porte), ...
            ecart_K, diff_porte, diff_adapte, duree, ternaire(ok, 'identique', 'DIFFERENT'));
    fprintf('       gains finaux : P %.6f, I %.4f, D %.6e (banc : %.6f, %.4f, %.6e)\n', K(end, 1), K(end, 2), ...
            K(end, 3), r.K(end, 1), r.K(end, 2), r.K(end, 3));
end
if tout_ok
    fprintf('\nRESULTAT : le bloc elm_pid_adaptatif.m reproduit banc_elm_pid.py. Etape suivante : Construction_ELM_PID.m.\n');
else
    error('Le bloc elm_pid_adaptatif.m ne reproduit pas le banc : ne pas construire le modele Simulink.');
end


%% ===================== Fonctions =====================================

function s = ternaire(condition, si_vrai, si_faux)
    % Choix entre deux textes.
    if condition
        s = si_vrai;
    else
        s = si_faux;
    end
end
