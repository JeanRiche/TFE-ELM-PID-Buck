function sc = charger_scenario(code)
% CHARGER_SCENARIO  Charge un essai commun dans le workspace de base.
%
% OBJECTIF
% --------
% Lire le fichier scenario_<code>.mat ecrit par scenarios_communs.py et
% creer, dans le workspace de base (celui ou Simulink cherche les
% variables des blocs), les variables qu'utilise le modele commun
% Buck_Commun.slx et tous les modeles construits a partir de lui :
%   sc_dvref  [t, dvref]  ajout a la consigne de 100 V (V)
%   sc_vin    [t, vin]    tension d'entree (V)
%   sc_gx     [t, gx]     conductance de la charge electronique (S)
%   sc_bruit  [t, bruit]  bruit ajoute a la mesure de Vout (V)
%   sc_R0     resistance physique de la charge (ohms)
%   sc_q      pas de quantification de la mesure (V ; 1e-9 = sans effet)
%   sc_duree  duree de l'essai (s), a passer a sim() comme StopTime
%   sc_code   code de l'essai (texte)
%
% POURQUOI LES INSTANTS SONT DECALES D'UN DEMI-PAS
% -------------------------------------------------
% Les blocs From Workspace du modele sont echantillonnes a la periode du
% regulateur Te = Tc = 1/220000 s, sans interpolation : a l'instant k*Te,
% ils sortent la derniere valeur dont l'instant est inferieur ou egal a
% k*Te. Si l'on donnait les instants
% exacts k*Te, un arrondi de 1e-20 s pourrait faire choisir la valeur du
% pas precedent au moment d'un saut. Avec les instants 0, Te/2, 3Te/2,
% 5Te/2, ..., la valeur numero k est choisie sans ambiguite a l'instant
% k*Te : c'est exactement ce que fait le banc Python au pas k.
%
% UTILISATION
% -----------
%   sc = charger_scenario('S9');
%   sim('Buck_Commun', 'StopTime', num2str(sc.duree));
% La structure renvoyee contient aussi le nom, la description, les
% instants des evenements et les profils bruts.
%
% Prerequis : scenario_<code>.mat dans le dossier courant MATLAB.
% Compatible MATLAB R2024a.
% ---------------------------------------------------------------------

    fichier = fullfile(pwd, ['scenario_' code '.mat']);
    if ~isfile(fichier)
        error('Fichier introuvable : %s (lancer d''abord scenarios_communs.py).', fichier);
    end
    m = load(fichier);                                  % variables ecrites par Python
    Te = double(m.Te);                                  % periode commune des profils (Tc = 1/220000 s)
    % Garde-fou : la periode des profils doit etre celle du regulateur.
    if ~(abs(Te - 1/220000) <= 1e-15)
        error(['%s : periode des profils %.6g s au lieu de 1/220000 s. Fichier de la ' ...
               'Utiliser les scenario_*.mat produits par scenarios_communs.py.'], ...
              ['scenario_' code '.mat'], Te);
    end
    n = numel(m.t);                                     % nombre de pas, instant 0 compris
    td = [0; ((1:n-1)' - 0.5) * Te];                    % instants decales d'un demi-pas
    assignin('base', 'sc_dvref', [td, double(m.dvref(:))]);
    assignin('base', 'sc_vin',   [td, double(m.vin(:))]);
    assignin('base', 'sc_gx',    [td, double(m.gx(:))]);
    assignin('base', 'sc_bruit', [td, double(m.bruit(:))]);
    assignin('base', 'sc_R0',    double(m.R0));
    assignin('base', 'sc_q',     double(m.q));
    assignin('base', 'sc_duree', double(m.duree));
    assignin('base', 'sc_code',  strtrim(char(m.code)));

    sc = struct();
    sc.code = strtrim(char(m.code));
    sc.nom = strtrim(char(m.nom));
    sc.description = strtrim(char(m.description));
    sc.duree = double(m.duree);
    sc.Te = Te;
    sc.R0 = double(m.R0);
    sc.q = double(m.q);
    sc.evenements = double(m.evenements(:))';           % ligne d'instants (peut etre vide)
    sc.t = double(m.t(:));
    sc.dvref = double(m.dvref(:));
    sc.vin = double(m.vin(:));
    sc.gx = double(m.gx(:));
    sc.bruit = double(m.bruit(:));
end
