classdef ordonnanceur_flou_zhao < matlab.System
    %% ORDONNANCEUR_FLOU_ZHAO -- ordonnanceur flou des gains du PID (Zhao et al., 1993)
    %
    % CE QUE FAIT CE BLOC
    % -------------------
    % A chaque periode Tc, il calcule les gains P, I, D du bloc "PID
    % Controller" (forme Parallel) a partir de l'erreur e(k) et de sa
    % difference De(k) = e(k) - e(k-1), par les regles floues de l'article
    % (tableaux I, II, III, eq. 6 a 11). Le PID lui-meme reste le bloc de
    % Ziegler-Nichols, dont P, I et D sont passes en entrees externes
    % (Construction_Fuzzy_PID.m). C'est la copie exacte de la fonction zhao
    % de banc_fuzzy_pid.py : memes nombres, meme ordre des operations.
    % Tester_Fuzzy_PID_Rejeu.m verifie que les deux donnent les memes gains
    % pas par pas.
    %
    % ENTREE ET SORTIE (une fois par periode Tc, imposee par un Zero-Order
    % Hold a 1/(22000*10) s place devant le bloc)
    %   entree 1 : e = consigne - mesure (sortie de Sum1, comme le PID)
    %   sortie 1 : G = [P I D] (1 x 3), vers les entrees P, I, D du PID
    %
    % L'ALGORITHME, A CHAQUE PERIODE Tc
    % ---------------------------------
    % 1. De = e - e(k-1) (0 au premier pas) ; e' = e / E_M, De' = De / DE_M.
    % 2. Appartenances aux sept ensembles NB NM NS ZO PS PM PB : triangles
    %    centres en -1, -2/3, ..., 1, demi-largeur 1/3, extremites
    %    prolongees.
    % 3. Pour chaque regle de poids non nul (ligne par ligne) :
    %    w = mu_e(i) mu_De(j) ; q = exp(-4 w) ;
    %    K'p += w (1 - q) si Big, w q si Small ; K'd de meme ; alpha += w alpha(i,j).
    % 4. P = KP_MIN + (KP_MAX - KP_MIN) K'p ; D = KD_MIN + (KD_MAX - KD_MIN) K'd ;
    %    I = P^2 / (alpha D).
    %
    % FICHIER LU (dans le dossier courant, un seul exemplaire)
    %   fuzzy_pid_reglages.mat : KP_MIN, KP_MAX, KD_MIN, KD_MAX, E_M, DE_M,
    %                            TABLE_KP, TABLE_KD (1 = Big, 0 = Small),
    %                            TABLE_ALPHA, K_ZN (banc_fuzzy_pid.py)
    % Aucun nombre de la methode n'est ecrit dans ce fichier.
    %
    % GainsFixes = true : sort les gains de Ziegler-Nichols a chaque pas.
    % Sert seulement a l'auto-test T0 de Construction_Fuzzy_PID.m (le modele
    % doit alors redonner le PID classique).
    %
    % Utilise dans Simulink en execution interpretee. Compatible MATLAB
    % R2024a. N'utilise pas la Fuzzy Logic Toolbox.

    properties (Nontunable)
        FichierReglages = 'fuzzy_pid_reglages.mat';  % reglages du bloc
        GainsFixes (1,1) logical = false;            % true : gains de Ziegler-Nichols (auto-test T0)
    end

    properties (Access = private)
        KPMIN; KPMAX; KDMIN; KDMAX; EM; DEM; TKP; TKD; TAL; KZN
        e1; premier; der_zh; npas
    end

    methods
        function obj = ordonnanceur_flou_zhao(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function zh = lire_zhao(obj)
            % [K'p K'd alpha] du dernier pas (pour le test de rejeu).
            zh = obj.der_zh;
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            if exist(obj.FichierReglages, 'file') ~= 2
                error('ordonnanceur_flou_zhao : %s introuvable (il doit etre dans le dossier du modele).', ...
                      obj.FichierReglages);
            end
            r = load(obj.FichierReglages);
            obj.KPMIN = double(r.KP_MIN); obj.KPMAX = double(r.KP_MAX);
            obj.KDMIN = double(r.KD_MIN); obj.KDMAX = double(r.KD_MAX);
            obj.EM = double(r.E_M); obj.DEM = double(r.DE_M);
            obj.TKP = double(r.TABLE_KP); obj.TKD = double(r.TABLE_KD); obj.TAL = double(r.TABLE_ALPHA);
            obj.KZN = double(r.K_ZN(:)');
            if ~isequal(size(obj.TKP), [7 7]) || ~isequal(size(obj.TKD), [7 7]) || ~isequal(size(obj.TAL), [7 7])
                error('ordonnanceur_flou_zhao : tables de taille inattendue dans %s.', obj.FichierReglages);
            end
            % Les deux regles citees dans le texte de l'article (points a1 et b1 de la fig. 4)
            if ~(obj.TKP(7, 4) == 1 && obj.TKD(7, 4) == 0 && obj.TAL(7, 4) == 2 && ...
                 obj.TKP(4, 1) == 0 && obj.TKD(4, 1) == 1 && obj.TAL(4, 1) == 5)
                error('ordonnanceur_flou_zhao : les tables ne redonnent pas les regles citees dans l''article.');
            end
        end

        function resetImpl(obj)
            obj.e1 = 0;
            obj.premier = true;
            obj.der_zh = [NaN NaN NaN];
            obj.npas = 0;
        end

        function G = stepImpl(obj, e)
            if obj.premier
                de = 0;                               % ecart 4 : De(0) = 0
                obj.premier = false;
            else
                de = e - obj.e1;
            end
            obj.e1 = e;
            obj.npas = obj.npas + 1;
            if obj.GainsFixes
                G = obj.KZN;
                obj.der_zh = [NaN NaN NaN];
                return
            end
            [kp, kd, al] = ordonnanceur_flou_zhao.zhao(e / obj.EM, de / obj.DEM, obj.TKP, obj.TKD, obj.TAL);
            P = obj.KPMIN + (obj.KPMAX - obj.KPMIN) * kp;
            D = obj.KDMIN + (obj.KDMAX - obj.KDMIN) * kd;
            I = P * P / (al * D);
            G = [P, I, D];
            obj.der_zh = [kp, kd, al];
        end

        function n = getNumInputsImpl(~)
            n = 1;
        end
        function n = getNumOutputsImpl(~)
            n = 1;
        end
        function n1 = getInputNamesImpl(~)
            n1 = 'e';
        end
        function n1 = getOutputNamesImpl(~)
            n1 = 'PID';
        end
        function s1 = getOutputSizeImpl(~)
            s1 = [1 3];
        end
        function d1 = getOutputDataTypeImpl(~)
            d1 = 'double';
        end
        function c1 = isOutputComplexImpl(~)
            c1 = false;
        end
        function f1 = isOutputFixedSizeImpl(~)
            f1 = true;
        end
    end

    methods (Static)

        function m = appartenances(x)
            % Degres d'appartenance de x (normalise) a NB NM NS ZO PS PM PB
            % (memes formules que appartenances dans le banc).
            xs = min(max(x, -1), 1);
            m = zeros(1, 7);
            for j = 1:7
                m(j) = max(0, 1 - abs(3 * xs - (j - 4)));
            end
        end

        function [kp, kd, al] = zhao(en, den, tkp, tkd, tal)
            % K'p, K'd et alpha pour e' = en et De' = den (eq. 8 et 10),
            % regles de poids nul sautees, parcours ligne par ligne.
            mE = ordonnanceur_flou_zhao.appartenances(en);
            mD = ordonnanceur_flou_zhao.appartenances(den);
            kp = 0; kd = 0; al = 0;
            for i = 1:7
                if mE(i) > 0
                    for j = 1:7
                        if mD(j) > 0
                            w = mE(i) * mD(j);
                            q = exp(-4 * w);
                            if tkp(i, j) == 1
                                kp = kp + w * (1 - q);
                            else
                                kp = kp + w * q;
                            end
                            if tkd(i, j) == 1
                                kd = kd + w * (1 - q);
                            else
                                kd = kd + w * q;
                            end
                            al = al + w * tal(i, j);
                        end
                    end
                end
            end
        end
    end
end
