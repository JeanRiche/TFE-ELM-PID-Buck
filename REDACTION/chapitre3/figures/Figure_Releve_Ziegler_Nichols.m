%% Figure_Releve_Ziegler_Nichols.m
% =========================================================================
%                                                                         %
%   RELEVE DE ZIEGLER-NICHOLS ET CALCUL DES GAINS DU PID CLASSIQUE         %
%   (figure 1 et sections 4, 5 et 7 du document "Reglage du regulateur     %
%    PID classique par la methode de Ziegler-Nichols")                     %
%                                                                         %
%   Convertisseur Buck du travail : L = 10 mH, C = 47 uF, R = 5 ohms,      %
%   Vin = 200 V, decoupage a 22 kHz.                                       %
%                                                                         %
% =========================================================================
%
% CE QUE PRODUIT CE SCRIPT
%   1. La reponse indicielle du convertisseur (modele moyen), calculee
%      de facon exacte.
%   2. Le releve de la methode de la reponse indicielle : point
%      d'inflexion, tangente, retard apparent La (note L chez Mudry),
%      instant t2 (63 %),
%      constante de temps apparente T, t3, pente normalisee p, a = p*La
%      et temps mort relatif tau.
%   3. Les gains du PID par la table 1 de la note de F. Mudry (ligne
%      PID) : Kp = 1,2/(p*La*K0), Ti = 2*La, Td = 0,5*La.
%   4. Les valeurs a saisir dans le bloc "PID Controller" de Simulink,
%      en forme Parallel (retenue) et en forme Ideal (pour comparaison),
%      avec le coefficient N du filtre derive et les controles qui le
%      justifient (pole du filtre sous fs/2, stabilite de la
%      discretisation Forward Euler).
%   5. La figure du releve, a inserer dans le memoire, en PNG et en PDF
%      vectoriel.
%
% NOTATION
%   Le retard apparent est note La (L dans la note de Mudry) pour ne pas
%   le confondre avec l'inductance L du convertisseur.
%
% POURQUOI UNE REPONSE CALCULEE PLUTOT QUE SIMULEE
%   Le modele moyen du convertisseur est un second ordre a deux poles
%   reels. Sa reponse indicielle a une expression exacte, et le point
%   d'inflexion aussi (annulation de la derivee seconde). Le calcul
%   exact ne depend ni d'un pas d'integration ni de la phase de la
%   porteuse MLI a l'instant de l'echelon. Le document compare ce
%   releve a celui de la reponse commutee (section 4.3) : La differe de
%   11,9 us (pres d'un quart de periode de decoupage) et Kp de 8,9 %
%   (tableau 6 ; K0 = 197,2 V pour la reponse commutee. Les 7,4 %
%   annonces auparavant gardaient K0 = 200 V pour cette reponse).
%
% SORTIES
%   figures_zn/fig3_02_releve_Ziegler_Nichols.pdf  (vectoriel, 12 cm de large,
%       nom de l'image du memoire : a deposer tel quel dans figures/ d'Overleaf)
%   figures_zn/fig3_02_releve_Ziegler_Nichols.png  (300 dpi, meme taille)
%   Mise en page : 12 cm x 6,5 cm (0,75 x 16 cm de largeur de texte),
%   Times New Roman (repli Helvetica), 9 pt (graduations, legende,
%   annotations), 10 pt (etiquettes d'axes), traits >= 0,75 pt, pas de titre.
%   Toutes les grandeurs sont affichees dans la fenetre de commandes.
%
% VALEURS ATTENDUES (a comparer a l'affichage)
%   La = 0.155951 ms, T = 1.857541 ms, p = 409.687 1/s, tau = 0.0775
%   Kp = 0.093910, Ti = 0.311901 ms, Td = 0.077975 ms
%   Parallel : P = 0.093910, I = 301.089, D = 7.3227e-06, N = 64122.9
%
% Prerequis : aucun (MATLAB de base, R2024a ou plus recent).
% =========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------
%  1. PARAMETRES DU CONVERTISSEUR ET DU REGULATEUR
%  ------------------------------------------------------------------
L_ind = 10e-3;           % inductance [H] (note de dimensionnement)
C_cap = 47e-6;           % capacite [F] (note de dimensionnement)
R_ch  = 5;               % charge nominale [ohm] : c'est le point de reglage
Vin   = 200;             % tension d'entree [V]
fs    = 22000;           % frequence de decoupage [Hz]

E     = 1;               % amplitude de l'echelon de rapport cyclique
                         % (sans unite). Le modele etant lineaire, la
                         % reponse normalisee y/y(inf) ne depend pas de E :
                         % un echelon de 0 a 0,5 donne le meme releve.

Tc    = 1/(fs*10);       % periode d'echantillonnage du regulateur [s] :
                         % dix echantillons par periode de decoupage,
                         % synchronises avec la porteuse (section 6)
N_mudry = 5;             % coefficient du filtre derive dans la convention
                         % de la note de Mudry (choix justifie en 7.2)

%% ------------------------------------------------------------------
%  2. MODELE MOYEN ET CONSTANTES DE TEMPS
%  ------------------------------------------------------------------
% Fonction de transfert du rapport cyclique vers la tension de sortie :
%       Gvd(s) = Vin / ( L*C*s^2 + (L/R)*s + 1 )
% roots() attend les coefficients par puissances DECROISSANTES.
den   = [L_ind*C_cap, L_ind/R_ch, 1];
poles = roots(den);                          % deux poles reels negatifs
tau   = sort(-1./real(poles), 'descend');    % constantes de temps [s]
tau1  = tau(1);                              % constante de temps lente
tau2  = tau(2);                              % constante de temps rapide

wn    = 1/sqrt(L_ind*C_cap);                 % pulsation propre [rad/s]
zeta  = (L_ind/R_ch)/(2*sqrt(L_ind*C_cap));  % amortissement (sans unite)

yinf  = Vin*E;                               % valeur finale [V]
K0    = yinf/E;                              % gain statique [V par unite]

%% ------------------------------------------------------------------
%  3. REPONSE INDICIELLE EXACTE ET POINT D'INFLEXION
%  ------------------------------------------------------------------
% Second ordre aperiodique a poles reels distincts, circuit au repos :
%   y(t)  = yinf*[1 - (tau1*exp(-t/tau1) - tau2*exp(-t/tau2))/(tau1 - tau2)]
%   y'(t) = yinf*(exp(-t/tau1) - exp(-t/tau2))/(tau1 - tau2)
% Le point d'inflexion annule y''(t), ce qui donne un instant exact :
%   t_infl = tau1*tau2*ln(tau1/tau2)/(tau1 - tau2)
y_fun  = @(t) yinf*(1 - (tau1*exp(-t/tau1) - tau2*exp(-t/tau2))/(tau1 - tau2));
dy_fun = @(t) yinf*(exp(-t/tau1) - exp(-t/tau2))/(tau1 - tau2);

t_infl = tau1*tau2*log(tau1/tau2)/(tau1 - tau2);   % instant d'inflexion [s]
y_infl = y_fun(t_infl);                            % tension a l'inflexion [V]
pente  = dy_fun(t_infl);                           % pente maximale [V/s]

% Controle independant : maximum de la derivee calculee numeriquement
% sur une grille de 1 ns. Le sommet de la derivee est tres plat : les deux
% instants coincident a quelques nanosecondes pres, sur 597 us.
dt   = 1e-9;
t    = 0:dt:6e-3;
y    = y_fun(t);
[~, idx] = max(gradient(y, dt));
ecart_ns = abs(t(idx) - t_infl)*1e9;

%% ------------------------------------------------------------------
%  4. GRANDEURS DU RELEVE
%  ------------------------------------------------------------------
% Tangente en (t_infl, y_infl) :  y = y_infl + pente*(t - t_infl).
L_ret = t_infl - y_infl/pente;                % t1 = La : coupe l'axe des temps
t3    = t_infl + (yinf - y_infl)/pente;       % t3 : coupe l'asymptote y(inf)

% t2 : instant ou la reponse atteint 63 % de sa valeur finale. On le
% resout exactement avec fzero, entre l'inflexion et 10*tau1.
t2_63 = fzero(@(x) y_fun(x) - 0.63*yinf, [t_infl, 10*tau1]);

T_app = t2_63 - L_ret;          % constante de temps apparente [s]
p_n   = pente/yinf;             % pente de la reponse NORMALISEE [1/s]
a_int = p_n*L_ret;              % ordonnee (normalisee) de la tangente en t = 0
tau_r = L_ret/(L_ret + T_app);  % temps mort relatif (sans unite)

%% ------------------------------------------------------------------
%  5. GAINS DE ZIEGLER-NICHOLS (table 1 de Mudry, ligne PID)
%  ------------------------------------------------------------------
Kp = 1.2/(p_n*L_ret*K0);        % gain proportionnel [unite de rapport cyclique par V]
Ti = 2*L_ret;                   % temps d'integration [s], Ti = 2*La
Td = 0.5*L_ret;                 % temps de derivation [s], Td = 0,5*La

% Champs du bloc PID de Simulink.
%   Forme Parallel : u = P*e + I*integrale(e) + D*N/(1 + N/s)*e
%   Forme Ideal    : u = P*[ e + I*integrale(e) + D*N/(1 + N/s)*e ]
P_par = Kp;      I_par = Kp/Ti;   D_par = Kp*Td;     % forme retenue
P_id  = Kp;      I_id  = 1/Ti;    D_id  = Td;        % pour comparaison

% Filtre du terme derive. Mudry : sTd/(1 + sTd/N), pole a N/Td.
% Simulink : pole a N (rad/s). Il faut donc saisir N_simulink = N_mudry/Td.
N_sim    = N_mudry/Td;                 % [rad/s]
f_filtre = N_sim/(2*pi);               % frequence du pole [Hz]
N_max    = pi*fs*Td;                   % plus grand N_mudry qui garde le pole sous fs/2

% Stabilite du filtre discretise en Forward Euler : son pole vaut
% z = 1 - N*Tc, stable si et seulement si 0 < N*Tc < 2.
NTc = N_sim*Tc;

%% ------------------------------------------------------------------
%  6. AFFICHAGE
%  ------------------------------------------------------------------
fprintf('=== Procede (modele moyen) ===\n');
fprintf('Poles                       : %.2f et %.2f 1/s\n', sort(real(poles),'descend'));
fprintf('Constantes de temps         : tau1 = %.4f ms , tau2 = %.5f ms\n', tau1*1e3, tau2*1e3);
fprintf('Pulsation propre            : %.1f rad/s (%.1f Hz)\n', wn, wn/(2*pi));
fprintf('Amortissement zeta          : %.3f (procede aperiodique)\n', zeta);
fprintf('Gain statique K0            : %.1f V par unite de rapport cyclique\n\n', K0);

fprintf('=== Releve sur la reponse indicielle ===\n');
fprintf('Point d''inflexion           : t = %.5f ms , y = %.3f V\n', t_infl*1e3, y_infl);
fprintf('Controle numerique          : ecart de %.1f ns sur l''instant d''inflexion\n', ecart_ns);
fprintf('Pente au point d''inflexion  : %.1f V/s (p normalisee = %.3f 1/s)\n', pente, p_n);
fprintf('Retard apparent La = t1     : %.6f ms\n', L_ret*1e3);
fprintf('t2 (63 %% de y(inf))         : %.5f ms\n', t2_63*1e3);
fprintf('Constante de temps T        : %.6f ms\n', T_app*1e3);
fprintf('t3                          : %.5f ms\n', t3*1e3);
fprintf('a = p*La                    : %.5f\n', a_int);
fprintf('Temps mort relatif tau      : %.4f\n\n', tau_r);

fprintf('=== Gains de Ziegler-Nichols (forme standard) ===\n');
fprintf('Kp = %.6f    Ti = %.6f ms    Td = %.6f ms\n\n', Kp, Ti*1e3, Td*1e3);

fprintf('=== Champs du bloc PID Controller ===\n');
fprintf('Forme Parallel (retenue) : P = %.6f  I = %.3f  D = %.4e\n', P_par, I_par, D_par);
fprintf('Forme Ideal              : P = %.6f  I = %.3f  D = %.4e\n', P_id, I_id, D_id);
fprintf('N (Simulink)             : %.1f rad/s, pole a %.0f Hz\n', N_sim, f_filtre);
fprintf('Controle fs/2            : N_mudry = %d, plus grand N admis = %.2f\n', N_mudry, N_max);
fprintf('Controle Forward Euler   : N*Tc = %.3f (stable si < 2)\n', NTc);
fprintf('Sample time              : 1/(22000*10) = %.4e s\n', Tc);
if N_mudry > N_max
    warning('Le pole du filtre derive est au-dessus de fs/2.');
end
if NTc >= 2
    warning('Filtre derive instable en Forward Euler avec ce pas.');
end

%% ------------------------------------------------------------------
%  7. FIGURE DU RELEVE
%  ------------------------------------------------------------------
% Taille d'insertion dans le memoire : 12 cm de large, proportions de
% l'ancienne figure (1000 x 540 px).
LARG_CM = 12;  HAUT_CM = LARG_CM * 540/1000;
POLICE  = 'Times New Roman';
if ~any(strcmpi(listfonts, POLICE)); POLICE = 'Helvetica'; end
fig = figure('Units', 'centimeters', 'Position', [5 5 LARG_CM HAUT_CM], 'Color', 'w');
hold on; grid on; box on;

k_aff = 1:100:numel(t);                                          % un point sur cent suffit au trace
h1 = plot(t(k_aff)*1e3, y(k_aff), 'b', 'LineWidth', 1.8);       % reponse
h2 = plot([L_ret t3]*1e3, [0 yinf], 'r--', 'LineWidth', 1.4);   % tangente, de (La,0) a (t3,yinf)

yline(yinf,      ':', 'Color', [0 0 0],       'LineWidth', 0.9);
yline(0.63*yinf, ':', 'Color', [.45 .45 .45], 'LineWidth', 0.9);

h3 = plot(t_infl*1e3, y_infl,    'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
h4 = plot(L_ret*1e3,  0,         'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
h5 = plot(t2_63*1e3,  0.63*yinf, 'kd', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
h6 = plot(t3*1e3,     yinf,      'k^', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

% Cote de la constante de temps apparente T, dessinee en coordonnees
% des donnees (elle reste en place si la fenetre est redimensionnee).
y_cote = -14;
plot([L_ret t2_63]*1e3, [y_cote y_cote], '-', 'Color', [0 .5 0], 'LineWidth', 1.4);
plot([L_ret L_ret]*1e3, y_cote + [-5 5],  '-', 'Color', [0 .5 0], 'LineWidth', 1.4);
plot([t2_63 t2_63]*1e3, y_cote + [-5 5],  '-', 'Color', [0 .5 0], 'LineWidth', 1.4);
text((L_ret + t2_63)*1e3/2, y_cote - 12, sprintf('T = t_2 - t_1 = %.4f ms', T_app*1e3), ...
     'Color', [0 .5 0], 'HorizontalAlignment', 'center', 'FontSize', 9);

% Etiquettes des niveaux a gauche, au-dessus des pointilles, la ou ni la
% reponse ni la tangente ne passent.
text(0.05, yinf + 9,      sprintf('y(\\infty) = %.0f V', yinf), ...
     'FontSize', 9, 'HorizontalAlignment', 'left');
% 63 % : a droite, juste au-dessus des pointilles (la reponse y passe
% vers 185 V, la tangente est deja terminee, la legende est en dessous).
text(5.95, 0.63*yinf + 9, sprintf('63 %% de y(\\infty) = %.0f V', 0.63*yinf), ...
     'FontSize', 9, 'Color', [.45 .45 .45], 'HorizontalAlignment', 'right');

xlabel('Temps (ms)', 'FontSize', 10);
ylabel('Tension de sortie V_{out} (V)', 'FontSize', 10);
% Titre retire : la legende de la figure est dans le memoire (figure 3.2).
% title({'Releve de Ziegler-Nichols sur la reponse indicielle en boucle ouverte', ...
%        sprintf('(modele moyen, L = 10 mH, C = 47 \\muF, R = 5 \\Omega, echelon E = %g)', E)});
legend([h1 h2 h3 h4 h5 h6], ...
       {'Réponse indicielle du convertisseur', ...
        'Tangente au point d''inflexion', ...
        'Point d''inflexion', ...
        sprintf('t_1 = L_a = %.4f ms', L_ret*1e3), ...
        sprintf('t_2 = %.4f ms  (63 %%)', t2_63*1e3), ...
        sprintf('t_3 = %.4f ms', t3*1e3)}, ...
       'Location', 'southeast', 'FontSize', 9, 'Box', 'off');
xlim([-0.1 6]);
ylim([-38 225]);

%% ------------------------------------------------------------------
%  8. ENREGISTREMENT (PDF vectoriel sous le nom du memoire + PNG 300 dpi)
%  ------------------------------------------------------------------
set(findall(fig, '-property', 'FontName'), 'FontName', POLICE);
set(gca, 'FontSize', 9, 'Layer', 'top');       % graduations 9 pt
set(get(gca, 'XLabel'), 'FontSize', 10);       % etiquettes d'axes 10 pt
set(get(gca, 'YLabel'), 'FontSize', 10);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [LARG_CM HAUT_CM], ...
    'PaperPosition', [0 0 LARG_CM HAUT_CM]);
drawnow;
% Dossier de sortie : figures_zn, a cote de ce script (quel que soit le
% dossier courant de MATLAB).
ici = fileparts(mfilename('fullpath'));
if isempty(ici); ici = pwd; end
dossier_sortie = fullfile(ici, 'figures_zn');
if ~exist(dossier_sortie, 'dir'); mkdir(dossier_sortie); end
base = fullfile(dossier_sortie, 'fig3_02_releve_Ziegler_Nichols');
try
    exportgraphics(fig, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
catch
    print(fig, [base '.png'], '-dpng', '-r300');
end
if exist('exportgraphics', 'file')
    try
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', ...
                       'BackgroundColor', 'white', 'Padding', 'figure');
    catch
        exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
else
    set(fig, 'InvertHardcopy', 'off');
    print(fig, [base '.pdf'], '-dpdf', '-vector');
end
if exist([base '.pdf'], 'file')
    fprintf('\nFigure enregistree :\n  %s.pdf\n  %s.png\n', base, base);
else
    warning('Le PDF n''a pas ete cree : %s.pdf', base);
end
