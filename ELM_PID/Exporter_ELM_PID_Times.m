%% Exporter_ELM_PID_Times.m
% =========================================================================
%                                                                         %
%   EXPORT PDF DU MODELE SIMULINK DE L'ELM-PID, EN TIMES NEW ROMAN         %
%   (figure 4.6 du memoire : fig4_06_modele_Simulink_ELM_PID.pdf)          %
%                                                                         %
% =========================================================================
%
% CE QUE FAIT CE SCRIPT
%   1. Ouvre Buck_Commun_ELM_PID.slx, situe dans le dossier de ce script.
%      S'il n'existe pas, le construit d'abord en lancant
%      Construction_ELM_PID.m (deux a quatre minutes, auto-tests compris),
%      sans modifier ce script. Si les auto-tests echouent, l'erreur est
%      remontee ; le .slx deja enregistre sera utilise a la relance.
%   2. Met en Times New Roman, la police du texte du memoire, la police de
%      tous les blocs, de toutes les lignes (etiquettes de signaux) et de
%      toutes les annotations du premier niveau du schema.
%   3. Essaie une disposition automatique du schema
%      (Simulink.BlockDiagram.arrangeSystem) pour limiter les fils qui
%      traversent les etiquettes. Drapeau DISPOSITION_AUTO ci-dessous :
%      false pour la desactiver. Fonction absente des releases anciennes de
%      Simulink : l'essai est alors ignore, sans erreur.
%   4. Exporte le premier niveau du schema en PDF vectoriel dans le dossier
%      du script : fig4_06_modele_Simulink_ELM_PID.pdf, a la taille du schema
%      (PaperPositionMode = auto). Le chemin complet du fichier est affiche.
%   5. Referme le modele SANS l'enregistrer : Buck_Commun_ELM_PID.slx reste
%      tel quel (police et disposition ne changent que le temps de l'export).
%
% PREREQUIS
%   MATLAB avec Simulink et Simscape (les blocs du modele) ; police
%   Times New Roman installee ; tous les fichiers du dossier ELM_PID (ceux
%   que demande Construction_ELM_PID.m, dont elm_pid_adaptatif.m, necessaire
%   pour afficher le bloc "ELM-PID Adaptatif"). Le modele ne doit pas etre
%   deja ouvert : le script le referme sans enregistrer, ce qui perdrait des
%   modifications en cours. Le dossier courant MATLAB est change pendant
%   l'execution, puis retabli.
%
% COMMENT LANCER
%   >> Exporter_ELM_PID_Times
% =========================================================================

% --- Reglages ------------------------------------------------------------
DISPOSITION_AUTO     = true;    % true : essaie arrangeSystem avant l'export ; false : schema tel quel
DISPOSITION_COMPLETE = false;   % true : reorganise tout le schema ('FullLayout') ; false : ajustement local

dossier = fileparts(mfilename('fullpath'));
modele  = 'Buck_Commun_ELM_PID';
fichier = fullfile(dossier, [modele '.slx']);
sortie  = fullfile(dossier, 'fig4_06_modele_Simulink_ELM_PID.pdf');
police  = 'Times New Roman';
taille  = 16;   % taille des textes (pt) : la figure est reduite d'environ moitie dans le memoire

if bdIsLoaded(modele)
    error('Exporter_ELM_PID_Times:ouvert', ...
          ['Le modele %s est deja ouvert. Le fermer (en enregistrant si besoin) ' ...
           'avant de lancer ce script.'], modele);
end

ancien = pwd;                       % dossier courant a retablir a la fin
cd(dossier);                        % Construction_ELM_PID.m travaille dans le dossier courant
try
    % Construction du modele s'il n'existe pas encore.
    if ~isfile(fichier)
        fprintf('%s.slx absent : construction par Construction_ELM_PID.m (deux a quatre minutes).\n', modele);
        construire_modele();
        if ~isfile(fichier)
            error('Exporter_ELM_PID_Times:construction', ...
                  'Construction_ELM_PID.m n''a pas produit %s.', fichier);
        end
    end

    open_system(fichier);

    % Blocs, lignes et annotations du premier niveau du schema.
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
    if ~isfile(sortie)
        error('Exporter_ELM_PID_Times:export', 'Le PDF n''a pas ete ecrit : %s', sortie);
    end
    fprintf('Ecrit : %s\n', sortie);
catch erreur
    if bdIsLoaded(modele)
        close_system(modele, 0);    % 0 : sans enregistrer
    end
    cd(ancien);
    rethrow(erreur);
end
close_system(modele, 0);            % 0 : sans enregistrer
cd(ancien);


%% ===================== Fonctions =====================================
% (a la fin du fichier, comme dans Construction_ELM_PID.m)

function construire_modele()
    % Lance Construction_ELM_PID.m dans l'espace de travail de CETTE
    % fonction : le "clear" du debut de ce script ne vide alors que cet
    % espace, et pas les variables de Exporter_ELM_PID_Times.
    Construction_ELM_PID;
end
