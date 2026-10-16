%% CONSTRUCTION_PSO_PID.m
%
% OBJECTIF
% --------
% Construire Buck_Commun_PSO_PID.slx, le modele du PSO-PID (essaim de Gaing 2004 ;
% reglage hors ligne). C'est une copie de Buck_Commun.slx dans laquelle on
% change seulement les gains P, I et D du bloc "PID Controller". Les gains
% sont lus dans predictions_banc_pso_pid.json, ecrit par
% banc_pso_pid.py : Simulink utilise donc exactement les gains que le
% banc a simules.
%
% CE QUE LE SCRIPT CHANGE, ET RIEN D'AUTRE
% -----------------------------------------
%   P, I et D du bloc "PID Controller", ecrits avec 17 chiffres
%   significatifs. N (64122.9), la periode 1/(22000*10), la forme Parallel,
%   le Forward Euler, les bornes [0.01 ; 0.99], le clamping et les
%   conditions initiales ne bougent pas. Buck_Commun.slx n'est jamais
%   modifie : on part toujours d'une copie fraiche.
%
% CONTROLES
% ---------
%   - avant : le bloc PID de Buck_Commun.slx a bien les gains de
%     Ziegler-Nichols et le N attendu (on part du bon modele) ;
%   - apres : les gains relus dans le fichier enregistre valent ceux du
%     fichier JSON a 1e-12 pres en relatif.
%
% ORDRE D'EXECUTION (PSO-PID)
% --------------------------
%   1. recherche_pso_pid.py ; 2. banc_pso_pid.py ; 3. ce script ;
%   4. Simuler_PSO_PID.m.
%
% Prerequis dans le dossier courant MATLAB : Buck_Commun.slx (construit et
% verifie) et predictions_banc_pso_pid.json.
% Duree : quelques secondes. Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

clear; clc;

SOURCE  = 'Buck_Commun';                          % modele commun (jamais modifie)
MDL     = 'Buck_Commun_PSO_PID';                 % modele construit
BLOC    = 'PID Controller';                       % bloc dont on change les gains
JSON    = 'predictions_banc_pso_pid.json';       % gains choisis par banc_pso_pid.py
ZN      = struct('P', 0.093910, 'I', 301.089, 'D', 7.3227e-6, 'N', 64122.9);   % gains de depart attendus

% --- Prerequis ---
for f = {[SOURCE '.slx'], JSON}
    if ~isfile(fullfile(pwd, f{1}))
        error('Fichier introuvable dans le dossier courant : %s', f{1});
    end
end
pred = jsondecode(fileread(fullfile(pwd, JSON)));
if ~isfield(pred, 'format') || pred.format ~= 2 || ~isfield(pred, 'gains')
    error('%s ne vient pas de banc_pso_pid.py.', JSON);
end
g = pred.gains;                                   % P, I, D, N du PSO-PID
fprintf('Gains lus dans %s : P = %.10g, I = %.10g, D = %.10g, N = %.10g\n', JSON, g.P, g.I, g.D, g.N);
if abs(g.N - ZN.N) > 1e-9 * ZN.N
    error('Le N du fichier JSON (%.10g) n''est pas celui du PID de reference (%.10g).', g.N, ZN.N);
end

% --- Modeles ouverts : on ne perd jamais une modification non enregistree ---
if bdIsLoaded(SOURCE) && strcmp(get_param(SOURCE, 'Dirty'), 'on')
    error('%s est ouvert avec des modifications non enregistrees : l''enregistrer ou le fermer.', SOURCE);
end
for m = {MDL, SOURCE}
    if bdIsLoaded(m{1})
        close_system(m{1}, 0);                    % fermer sans enregistrer
    end
end

% --- Copie fraiche ---
chemin = fullfile(pwd, [MDL '.slx']);
if isfile(chemin)
    delete(chemin);
    fprintf('Ancien %s.slx supprime (on repart toujours d''une copie fraiche).\n', MDL);
end
copyfile(fullfile(pwd, [SOURCE '.slx']), chemin);
copies = which([MDL '.slx'], '-all');
if numel(copies) ~= 1
    error('Il doit exister exactement un %s.slx sur le chemin MATLAB, il y en a %d.', MDL, numel(copies));
end
load_system(chemin);
bloc = [MDL '/' BLOC];
if getSimulinkBlockHandle(bloc) == -1
    close_system(MDL, 0);
    error('Bloc introuvable : %s', bloc);
end

% --- Controle avant : on part bien du PID de Ziegler-Nichols ---
noms = {'P', 'I', 'D', 'N'};
for i = 1:numel(noms)
    v = evalin('base', get_param(bloc, noms{i}));
    if ~(abs(v - ZN.(noms{i})) <= 1e-9 * abs(ZN.(noms{i})))
        close_system(MDL, 0);
        error('%s de %s vaut %.10g au lieu de %.10g : Buck_Commun.slx n''a pas le PID de reference.', ...
              noms{i}, SOURCE, v, ZN.(noms{i}));
    end
end

% --- Changement des trois gains ---
set_param(bloc, 'P', sprintf('%.17g', g.P), 'I', sprintf('%.17g', g.I), 'D', sprintf('%.17g', g.D));
save_system(MDL);
close_system(MDL, 0);

% --- Controle apres : gains relus dans le fichier enregistre ---
load_system(chemin);
for i = 1:3
    v = evalin('base', get_param(bloc, noms{i}));
    if ~(abs(v - g.(noms{i})) <= 1e-12 * abs(g.(noms{i})))
        close_system(MDL, 0);
        error('%s relu dans %s.slx vaut %.17g au lieu de %.17g.', noms{i}, MDL, v, g.(noms{i}));
    end
end
fprintf('Gains relus dans %s.slx : P = %s, I = %s, D = %s, N = %s\n', MDL, get_param(bloc, 'P'), ...
        get_param(bloc, 'I'), get_param(bloc, 'D'), get_param(bloc, 'N'));
close_system(MDL, 0);

fprintf(['\n%s.slx est construit. Lancer ensuite ' ...
         'Simuler_PSO_PID.m.\n'], MDL);
