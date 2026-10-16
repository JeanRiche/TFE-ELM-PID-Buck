%% Exporter_Buck_Commun_Times.m
% =========================================================================
%                                                                         %
%   EXPORT PDF DU MODELE SIMULINK COMMUN, EN TIMES NEW ROMAN               %
%   (figure 3.3 du memoire : fig3_03_modele_Simulink_commun.pdf)           %
%                                                                         %
% =========================================================================
%
% CE QUE FAIT CE SCRIPT
%   1. Ouvre Buck_Commun.slx, situe dans le dossier de ce script.
%   2. Met en Times New Roman, la police du texte du memoire, la police de
%      tous les blocs, de toutes les lignes (etiquettes de signaux) et de
%      toutes les annotations du schema.
%   3. Essaie une disposition automatique du schema
%      (Simulink.BlockDiagram.arrangeSystem) pour limiter les fils qui
%      traversent les etiquettes. Drapeau DISPOSITION_AUTO ci-dessous :
%      false pour la desactiver. Fonction absente des releases anciennes de
%      Simulink : l'essai est alors ignore, sans erreur.
%   4. Exporte le schema en PDF vectoriel dans le dossier du script :
%      Buck_Commun.pdf, a la taille du schema (PaperPositionMode = auto).
%   5. Referme le modele SANS l'enregistrer : Buck_Commun.slx reste tel quel
%      (police et disposition ne changent que le temps de l'export).
%
% FIGURE DU MEMOIRE
%   Buck_Commun.pdf est la figure 3.3 du memoire
%   (fig3_03_modele_Simulink_commun.pdf).
%
% PREREQUIS
%   MATLAB avec Simulink et Simscape (les blocs du modele) ; police
%   Times New Roman installee. Le modele ne doit pas etre deja ouvert : le
%   script le referme sans enregistrer, ce qui perdrait des modifications
%   en cours.
%
% COMMENT LANCER
%   >> Exporter_Buck_Commun_Times
% =========================================================================

% --- Reglages ------------------------------------------------------------
DISPOSITION_AUTO     = true;    % true : essaie arrangeSystem avant l'export ; false : schema tel quel
DISPOSITION_COMPLETE = false;   % true : reorganise tout le schema ('FullLayout') ; false : ajustement local

dossier = fileparts(mfilename('fullpath'));
modele  = 'Buck_Commun';
fichier = fullfile(dossier, [modele '.slx']);
sortie  = fullfile(dossier, [modele '.pdf']);
police  = 'Times New Roman';
taille  = 16;   % taille des textes (pt) : la figure est reduite d'environ moitie dans le memoire

if ~isfile(fichier)
    error('Exporter_Buck_Commun_Times:modele', 'Modele introuvable : %s', fichier);
end
if bdIsLoaded(modele)
    error('Exporter_Buck_Commun_Times:ouvert', ...
          ['Le modele %s est deja ouvert. Le fermer (en enregistrant si besoin) ' ...
           'avant de lancer ce script.'], modele);
end

open_system(fichier);
try
    % Blocs, lignes et annotations du niveau principal du schema.
    types = {'block', 'line', 'annotation'};
    noms  = {'blocs', 'lignes', 'annotations'};
    for k = 1:numel(types)
        h = find_system(modele, 'FindAll', 'on', 'SearchDepth', 1, 'Type', types{k});
        modifies = 0;
        for i = 1:numel(h)
            try
                set_param(h(i), 'FontName', police);
                try, set_param(h(i), 'FontSize', taille); catch, end
                modifies = modifies + 1;
            catch
                % objet sans propriete FontName : laisse tel quel
            end
        end
        fprintf('%-12s : %d trouves, %d en %s\n', noms{k}, numel(h), modifies, police);
    end

    % Disposition automatique (facultative, ignoree si la fonction est absente).
    if DISPOSITION_AUTO
        try
            if DISPOSITION_COMPLETE
                Simulink.BlockDiagram.arrangeSystem(modele, 'FullLayout', 'true');
            else
                Simulink.BlockDiagram.arrangeSystem(modele);
            end
            fprintf('Disposition automatique : appliquee.\n');
        catch avis
            fprintf('Disposition automatique : ignoree (%s).\n', avis.message);
        end
    end

    % PDF vectoriel, page a la taille du schema.
    set_param(modele, 'PaperPositionMode', 'auto');
    print(['-s' modele], '-dpdf', sortie);
    fprintf('Ecrit : %s\n', sortie);
catch erreur
    close_system(modele, 0);        % 0 : sans enregistrer
    rethrow(erreur);
end
close_system(modele, 0);            % 0 : sans enregistrer
