%% TESTER_FUZZY_PID_REJEU.m
%
% OBJECTIF
% --------
% Verifier, avant toute simulation Simulink, que le bloc MATLAB
% ordonnanceur_flou_zhao.m calcule exactement les gains que calcule
% banc_fuzzy_pid.py. Le banc a enregistre, pas par pas, l'erreur e recue
% par le regulateur sur quatre essais (S1, S7a, S8b, S9) et les gains
% [P I D] produits, ainsi que le systeme flou seul sur une grille de 255
% points (e, De). Ce script donne les memes entrees au bloc, hors
% Simulink, et compare. Le PID lui-meme est le bloc Simulink de
% Ziegler-Nichols : il est controle par l'auto-test T0 de
% Construction_Fuzzy_PID.m, pas ici.
%
% CE QUI EST COMPARE
% ------------------
%   - le systeme flou seul (K'p, K'd, alpha, puis P, I, D) sur la grille ;
%   - pour chaque essai et chaque pas : P, I, D et K'p, K'd, alpha.
% Tolerance : 1e-12 en relatif (les ecarts attendus viennent seulement des
% arrondis de exp, quelques 1e-16).
%
% ORDRE D'EXECUTION
% -----------------
%   1. banc_fuzzy_pid.py ; 2. ce script ; 3. Construction_Fuzzy_PID.m ;
%   4. Simuler_Fuzzy_PID.m.
%
% Prerequis dans le dossier courant MATLAB : ordonnanceur_flou_zhao.m,
% fuzzy_pid_reglages.mat, reference_rejeu_fuzzy_pid.mat.
% Duree : une demi-minute a deux minutes. Compatible R2024a.
% N'utilise pas la Fuzzy Logic Toolbox.
% ---------------------------------------------------------------------

clear; clc;

% --- Prerequis ---
for f = {'ordonnanceur_flou_zhao.m', 'fuzzy_pid_reglages.mat', 'reference_rejeu_fuzzy_pid.mat'}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
    if numel(which(f{1}, '-all')) ~= 1
        error('Il doit exister exactement un %s sur le chemin MATLAB.', f{1});
    end
end
clear ordonnanceur_flou_zhao                       % recharger la classe si elle a change

reg = load('fuzzy_pid_reglages.mat');
ref = load('reference_rejeu_fuzzy_pid.mat');

% --- 1. Le systeme flou seul ---
pts = double(ref.FIS_ENTREES);
zh_att = double(ref.FIS_SORTIES);
g_att = double(ref.FIS_GAINS);
ecart_zh = 0; ecart_g = 0;
for i = 1:size(pts, 1)
    [kp, kd, al] = ordonnanceur_flou_zhao.zhao(pts(i, 1) / reg.E_M, pts(i, 2) / reg.DE_M, ...
        double(reg.TABLE_KP), double(reg.TABLE_KD), double(reg.TABLE_ALPHA));
    P = reg.KP_MIN + (reg.KP_MAX - reg.KP_MIN) * kp;
    D = reg.KD_MIN + (reg.KD_MAX - reg.KD_MIN) * kd;
    I = P * P / (al * D);
    ecart_zh = max(ecart_zh, max(abs([kp kd al] - zh_att(i, :)) ./ max(abs(zh_att(i, :)), 1e-300)));
    ecart_g = max(ecart_g, max(abs([P I D] - g_att(i, :)) ./ abs(g_att(i, :))));
end
fprintf('Systeme flou seul : %d points, ecart relatif maximal %.1e sur (K''p, K''d, alpha), %.1e sur (P, I, D).\n', ...
        size(pts, 1), ecart_zh, ecart_g);
if ecart_zh > 1e-12 || ecart_g > 1e-12
    error('Le systeme flou du bloc ne redonne pas celui de banc_fuzzy_pid.py.');
end

% --- 2. Rejeu pas par pas ---
codes = {'S1', 'S7a', 'S8b', 'S9'};
tout_ok = true;
for c = 1:numel(codes)
    r = ref.(codes{c});
    e = double(r.e(:));
    N = numel(e);
    bloc = ordonnanceur_flou_zhao();
    G = zeros(N, 3);
    ZH = zeros(N, 3);
    t0 = tic;
    for k = 1:N
        G(k, :) = bloc(e(k));
        ZH(k, :) = lire_zhao(bloc);
    end
    duree = toc(t0);
    release(bloc);
    ecart_G = max(max(abs(G - double(r.K)) ./ abs(double(r.K))));
    ecart_Z = max(max(abs(ZH - double(r.ZH)) ./ max(abs(double(r.ZH)), 1e-300)));
    ok = ecart_G <= 1e-12 && ecart_Z <= 1e-12 && size(r.K, 1) == N;
    tout_ok = tout_ok && ok;
    fprintf('%-4s : %6d pas ; ecart relatif maximal sur P, I, D %.1e, sur K''p, K''d, alpha %.1e ; %4.1f s  -> %s\n', ...
            codes{c}, N, ecart_G, ecart_Z, duree, ternaire(ok, 'identique', 'DIFFERENT'));
    fprintf('       gains finaux : P %.6f, I %.4f, D %.4e (banc : %.6f, %.4f, %.4e)\n', G(end, 1), G(end, 2), ...
            G(end, 3), r.K(end, 1), r.K(end, 2), r.K(end, 3));
end
if tout_ok
    fprintf('\nRESULTAT : le bloc ordonnanceur_flou_zhao.m reproduit banc_fuzzy_pid.py. Etape suivante : Construction_Fuzzy_PID.m.\n');
else
    error('Le bloc ordonnanceur_flou_zhao.m ne reproduit pas le banc : ne pas construire le modele Simulink.');
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
