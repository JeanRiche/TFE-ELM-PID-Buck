classdef pinn_pid_adaptatif < matlab.System
    %% PINN_PID_ADAPTATIF -- adaptateur PINN des gains du PID
    %
    % Le bloc ne calcule pas la commande. Il sort les gains [P I D] du bloc
    % "PID Controller" de la base commune, passe en gains externes (meme
    % montage que l'ELM-PID et le Fuzzy-PID). Loi, cout et optimisation
    % d'Ito et Wasa avec la correction C1 (optimisation seulement si
    % l'erreur efficace de la fenetre depasse le seuil).
    %
    % CE QUE FAIT CE BLOC
    % -------------------
    % Toutes les 0.5 ms, si l'erreur efficace de la fenetre depasse SEUIL, il
    % reoptimise les multiplicateurs x des gains de Ziegler-Nichols
    % (K = K_DEPART .* x) par N_ADAM iterations d'Adam, depart a chaud. Le
    % gradient vient d'une simulation de NH pas (horizon) avec le bloc PID
    % parallele et un PINN qui predit le convertisseur, a partir de l'etat
    % estime par un observateur (filtre de Kalman etendu a charge nominale)
    % et des etats du bloc PID. C'est la copie exacte de la classe
    % AdaptateurPINN de banc_pinn_pid.py : memes reglages, memes noms, meme
    % ordre des operations. Tester_PINN_PID_Rejeu.m verifie que les deux
    % donnent les memes gains pas par pas.
    %
    % ENTREES ET SORTIE (une fois par periode Tc, imposee par trois
    % Zero-Order Hold a 1/(22000*10) s places devant le bloc)
    %   entree 1 : e      = consigne - mesure (sortie de Sum1)
    %   entree 2 : mesure = tension mesuree (jamais 100 - e)
    %   entree 3 : u      = sortie du bloc PID (rapport cyclique borne)
    %   sortie 1 : K      = gains [P I D] pour ce pas (1 x 3), vers le PID
    % La sortie ne depend que de l'etat du bloc (aucune traversee directe :
    % outputImpl / updateImpl). Au pas k, le bloc sort les gains, le PID
    % calcule u(k) avec eux, puis le bloc lit e(k), mesure(k), u(k). Il n'y
    % a donc pas de boucle algebrique entre le PID et ce bloc.
    %
    % L'ALGORITHME, A CHAQUE PERIODE Tc (mise a jour)
    % ------------------------------------------------
    % 0. Copie du bloc PID : avec e et les gains de ce pas, memes operations
    %    que le bloc (Forward Euler, filtre N, saturation [D_MIN ; D_MAX],
    %    clamping ; depart CI_INTEGRATEUR et CI_FILTRE) ; on en tire les
    %    etats xI (integrateur) et xF (filtre) apres le pas, dont part
    %    l'horizon. La commande recalculee est comparee a l'entree u (ecart
    %    affiche dans le bilan : il doit etre nul, ou de l'ordre de 1e-15).
    % 1. Observateur, correction par la mesure : etat [Vout iL Vin_eff],
    %    charge nominale G_NOM.
    % 2. Observateur, prediction d'un pas avec la fraction de conduction
    %    appliquee : s = min(max(ceil(u x NPER) - phase, 0), NREG) / NREG,
    %    phase = mod(k x NREG, NPER).
    % 3. Fenetre de NF periodes (0.5 ms) : somme des e^2.
    % 4. Fin de fenetre : si l'erreur efficace depasse SEUIL (correction C1),
    %    N_ADAM iterations d'Adam sur x depuis x courant, moments remis a
    %    zero ; cout de l'horizon
    %      J = moyenne de 1/2 (e^2 + W_U u^2) + W_F |x - CENTRE|^2,
    %    gradient par retropropagation dans le temps ; projection sur la
    %    boite [X_MIN ; X_MAX] apres chaque iteration ; K = K_DEPART .* x.
    %
    % FICHIERS LUS (dans le dossier courant, un seul exemplaire)
    %   pinn_pid_modele.mat   : poids du PINN (entrainement_pinn.py)
    %   pinn_pid_reglages.mat : boite, seuil, Adam, cout, observateur, bloc
    %                           PID, table du modele physique
    %                           (banc_pinn_pid.py)
    % Aucun reglage de l'algorithme n'est ecrit dans ce fichier ; seules
    % les normalisations d'entree du PINN (100, 60, 12, 10, 195, 45, comme
    % pinn_pas du banc) le sont, et les vecteurs de test les controlent.
    %
    % GainsFixes = true : le bloc sort les gains de Ziegler-Nichols a chaque
    % pas (pas d'optimisation) ; le modele doit alors redonner le PID
    % classique. Sert seulement aux auto-tests T0 des scripts de
    % construction.
    %
    % Utilise dans Simulink en execution interpretee (Construction_PINN_PID.m).
    % Compatible MATLAB R2024a.

    properties (Nontunable)
        FichierModele = 'pinn_pid_modele.mat';      % poids du PINN
        FichierReglages = 'pinn_pid_reglages.mat';  % reglages du bloc
        AfficherBilan (1,1) logical = true;         % une ligne de bilan a la fin de chaque simulation
        GainsFixes (1,1) logical = false;           % true : gains de Ziegler-Nichols (auto-tests T0)
    end

    properties (Access = private)
        % PINN (couches cachees tanh, puis couche lineaire)
        Wc; bc; nc; SVP; SIP; GC; GE
        % Reglages
        K0; XMIN; XMAX; CENTRE; SEUIL; NH; NADAM; ALPHA; BETA1; BETA2; EPSA; WU; WF; PSAT
        GNOM; VIN0; SIGY; SIGV; SIGI; SIGVIN; P0I; P0VIN
        NF; TC; NFILTRE; DMIN; DMAX; CI_I; CI_F
        NPER; NREG; LMOD; VF; PHI; GAM; NS; PRED_PINN
        % Gains et copie du bloc PID
        K; x; xI; xF; ecart_u
        % Observateur
        xo; Po; Qo; Ro; obs_init
        % Compteurs et fenetre
        k; cpt; se2; r; nfen; nb_adapt; der_adapte; der_J
    end

    methods
        function obj = pinn_pid_adaptatif(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function info = lire_fenetre(obj)
            % Etat de la derniere fenetre terminee (pour le test de rejeu).
            info = struct('nfen', obj.nfen, 'adapte', obj.der_adapte, 'J', obj.der_J, ...
                          'vin_eff', obj.xo(3), 'iL', obj.xo(2), 'ecart_u', obj.ecart_u);
        end

        function [vn, inn, D] = tester_pas_pinn(obj, v, i, s, vin)
            % Un pas du PINN et ses derivees (pour les vecteurs de test).
            [vn, inn, P, dvs, dis] = pas_pinn(obj, v, i, s, vin);
            D = [P(1,1) - 1, P(2,1); P(1,2), P(2,2) - 1; dvs, dis];
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            m = pinn_pid_adaptatif.charger(obj.FichierModele);
            r = pinn_pid_adaptatif.charger(obj.FichierReglages);
            if ~isfield(r, 'FORMAT') || double(r.FORMAT) ~= 3
                error(['pinn_pid_adaptatif : %s n''a pas le format ecrit par banc_pinn_pid.py. ' ...
                       'Ne pas melanger avec des fichiers d''une autre origine.'], ...
                      obj.FichierReglages);
            end
            if double(r.ENSEMBLE) ~= 0
                error('pinn_pid_adaptatif : %s demande la projection sur l''ensemble admissible (C3), absente de ce bloc.', ...
                      obj.FichierReglages);
            end
            % PINN : W1..Wn, b1..bn
            obj.Wc = {}; obj.bc = {};
            l = 1;
            while isfield(m, sprintf('W%d', l))
                obj.Wc{l} = double(m.(sprintf('W%d', l)));
                obj.bc{l} = double(m.(sprintf('b%d', l)));
                obj.bc{l} = obj.bc{l}(:)';
                l = l + 1;
            end
            obj.nc = numel(obj.Wc);
            if obj.nc < 2 || size(obj.Wc{1}, 1) ~= 6 || size(obj.Wc{end}, 2) ~= 2
                error('pinn_pid_adaptatif : tailles inattendues dans %s.', obj.FichierModele);
            end
            obj.SVP = double(m.SV); obj.SIP = double(m.SI);
            obj.GC = double(m.G_CENTRE); obj.GE = double(m.G_ECHELLE);
            if abs(double(m.TC) - double(r.TC)) > 1e-18
                error('pinn_pid_adaptatif : %s et %s n''ont pas la meme periode Tc.', obj.FichierModele, obj.FichierReglages);
            end
            % Reglages
            obj.K0 = double(r.K_DEPART(:)');
            obj.XMIN = double(r.X_MIN(:)'); obj.XMAX = double(r.X_MAX(:)');
            obj.CENTRE = double(r.CENTRE(:)');
            obj.SEUIL = double(r.SEUIL);                 % NaN : pas de seuil (optimisation a chaque fenetre)
            obj.NH = double(r.NH); obj.NADAM = double(r.N_ADAM);
            obj.ALPHA = double(r.ALPHA); obj.BETA1 = double(r.BETA1); obj.BETA2 = double(r.BETA2);
            obj.EPSA = double(r.EPS_ADAM);
            obj.WU = double(r.W_U); obj.WF = double(r.W_F); obj.PSAT = double(r.PENTE_SAT);
            obj.GNOM = double(r.G_NOM); obj.VIN0 = double(r.VIN_DEPART);
            obj.SIGY = double(r.SIG_Y); obj.SIGV = double(r.SIG_V); obj.SIGI = double(r.SIG_I);
            obj.SIGVIN = double(r.SIG_VIN); obj.P0I = double(r.P0_I); obj.P0VIN = double(r.P0_VIN);
            obj.NF = double(r.NF); obj.TC = double(r.TC); obj.NFILTRE = double(r.N_FILTRE);
            obj.DMIN = double(r.D_MIN); obj.DMAX = double(r.D_MAX);
            obj.CI_I = double(r.CI_INTEGRATEUR); obj.CI_F = double(r.CI_FILTRE);
            obj.NPER = double(r.NPER); obj.NREG = double(r.NREG);
            obj.LMOD = double(r.L_MOD); obj.VF = double(r.VF);
            obj.PHI = double(r.PHI); obj.GAM = double(r.GAM);          % (NS+1) x 4 ([p11 p12 p21 p22]) et (NS+1) x 2
            obj.NS = size(obj.PHI, 1) - 1;
            obj.PRED_PINN = double(r.PREDICTEUR_PINN) ~= 0;
            if numel(obj.K0) ~= 3 || numel(obj.XMIN) ~= 3 || numel(obj.XMAX) ~= 3 || numel(obj.CENTRE) ~= 3 ...
                    || size(obj.PHI, 2) ~= 4 || size(obj.GAM, 2) ~= 2 || size(obj.GAM, 1) ~= obj.NS + 1
                error('pinn_pid_adaptatif : tailles inattendues dans %s.', obj.FichierReglages);
            end
        end

        function resetImpl(obj)
            obj.x = [1 1 1];                          % multiplicateurs de Ziegler-Nichols
            obj.K = obj.K0 .* obj.x;
            obj.xI = obj.CI_I; obj.xF = obj.CI_F;     % copie des etats du bloc PID
            obj.ecart_u = 0;
            obj.obs_init = false;
            obj.xo = [0, 0, obj.VIN0]; obj.Po = zeros(3); obj.Qo = zeros(3); obj.Ro = 0;
            obj.k = 0; obj.cpt = 0; obj.se2 = 0; obj.r = 0;
            obj.nfen = 0; obj.nb_adapt = 0; obj.der_adapte = false; obj.der_J = NaN;
        end

        function K = outputImpl(obj, ~, ~, ~)
            % Gains du pas, d'apres l'etat seulement (pas de traversee directe).
            K = obj.K;
        end

        function updateImpl(obj, e, mesure, u)
            kk = obj.k;
            % 0. Copie du bloc PID avec les gains de ce pas (memes operations que BlocPID.pas du banc)
            Pg = obj.K(1); Ig = obj.K(2); Dg = obj.K(3);
            derivee = obj.NFILTRE * (Dg * e - obj.xF);
            b = Pg * e + obj.xI + derivee;
            ub = min(max(b, obj.DMIN), obj.DMAX);
            entree_I = Ig * e;
            if ~((b ~= ub) && (sign(entree_I) == sign(b - ub)))   % clamping
                obj.xI = obj.xI + obj.TC * entree_I;
            end
            obj.xF = obj.xF + obj.TC * derivee;
            obj.ecart_u = max(obj.ecart_u, abs(ub - u));
            if ~obj.obs_init                          % observateur : depart sur la premiere mesure
                obj.xo = [mesure, 0, obj.VIN0];
                if obj.SIGVIN > 0
                    p3 = obj.P0VIN;
                else
                    p3 = 0;
                end
                obj.Po = diag([obj.SIGY * obj.SIGY, obj.P0I, p3]);
                obj.Qo = diag([obj.SIGV * obj.SIGV, obj.SIGI * obj.SIGI, obj.SIGVIN * obj.SIGVIN]);
                obj.Ro = obj.SIGY * obj.SIGY;
                obj.obs_init = true;
            end
            % 1. Correction de l'observateur par la mesure
            P = obj.Po;
            S = P(1,1) + obj.Ro;
            Kk = P(:,1)' / S;                         % gain de Kalman (1 x 3)
            obj.xo = obj.xo + Kk * (mesure - obj.xo(1));
            obj.Po = P - Kk' * P(1,:);
            obj.r = e + mesure;                       % consigne courante
            % 2. Prediction de l'observateur avec la fraction de conduction appliquee
            ph = mod(kk * obj.NREG, obj.NPER);        % phase de la porteuse
            s_app = min(max(ceil(u * obj.NPER) - ph, 0), obj.NREG) / obj.NREG;
            [vn, inn, Pm, ~, ~, dv_vin, di_vin] = pas_modele(obj, obj.xo(1), obj.xo(2), s_app, obj.xo(3));
            F = [Pm(1,1), Pm(1,2), dv_vin; Pm(2,1), Pm(2,2), di_vin; 0, 0, 1];
            obj.xo = [vn, inn, obj.xo(3)];
            obj.Po = F * obj.Po * F' + obj.Qo;
            % 3. Fenetre
            obj.se2 = obj.se2 + e * e;
            obj.cpt = obj.cpt + 1;
            obj.k = obj.k + 1;
            if obj.cpt >= obj.NF
                eff = sqrt(obj.se2 / obj.NF);
                adapte = false;
                Jfin = NaN;
                if ~obj.GainsFixes && (isnan(obj.SEUIL) || eff > obj.SEUIL)   % C1 : zone morte
                    Jfin = optimiser(obj, obj.r);
                    adapte = true;
                    obj.nb_adapt = obj.nb_adapt + 1;
                end
                obj.der_adapte = adapte;
                obj.der_J = Jfin;
                obj.nfen = obj.nfen + 1;
                obj.cpt = 0; obj.se2 = 0;
            end
        end

        function [f1, f2, f3] = isInputDirectFeedthroughImpl(~, ~, ~, ~)
            f1 = false;
            f2 = false;
            f3 = false;
        end

        function releaseImpl(obj)
            if obj.AfficherBilan && ~isempty(obj.nfen) && obj.nfen > 0
                fprintf(['  PINN-PID : %d fenetres, optimisees sur %d ; gains finaux P=%.5f I=%.4f D=%.4e ' ...
                         '(x ZN %.3f %.3f %.3f) ; Vin_eff=%.1f V ; ecart maximal copie du PID / entree u %.1e\n'], ...
                        obj.nfen, obj.nb_adapt, obj.K(1), obj.K(2), obj.K(3), obj.x(1), obj.x(2), obj.x(3), ...
                        obj.xo(3), obj.ecart_u);
            end
        end

        function n = getNumInputsImpl(~)
            n = 3;
        end
        function n = getNumOutputsImpl(~)
            n = 1;
        end
        function [n1, n2, n3] = getInputNamesImpl(~)
            n1 = 'e';
            n2 = 'mesure';
            n3 = 'u';
        end
        function n1 = getOutputNamesImpl(~)
            n1 = 'K';
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

    methods (Access = private)

        function [vn, inn, P, dvs, dis, dv_vin, di_vin] = pas_modele(obj, v, i, s, vin)
            % Un pas Tc du modele physique nominal (table interpolee en s),
            % memes formules que modele_pas dans banc_pinn_pid.py.
            x = s * obj.NS;
            j = min(floor(x), obj.NS - 1); f = x - j;   % segment de la table (indice j+1 en MATLAB)
            a = obj.PHI(j + 1, :); b = obj.PHI(j + 2, :);
            Pl = a + f * (b - a); dPl = (b - a) * obj.NS;
            ga = obj.GAM(j + 1, :); gb = obj.GAM(j + 2, :);
            gg = ga + f * (gb - ga); dg = (gb - ga) * obj.NS;
            src = (s * vin - (1 - s) * obj.VF) / obj.LMOD; dsrc = (vin + obj.VF) / obj.LMOD;
            vn = Pl(1) * v + Pl(2) * i + gg(1) * src;
            inn = Pl(3) * v + Pl(4) * i + gg(2) * src;
            dvs = dPl(1) * v + dPl(2) * i + dg(1) * src + gg(1) * dsrc;
            dis = dPl(3) * v + dPl(4) * i + dg(2) * src + gg(2) * dsrc;
            P = [Pl(1), Pl(2); Pl(3), Pl(4)];
            dv_vin = gg(1) * s / obj.LMOD; di_vin = gg(2) * s / obj.LMOD;
        end

        function [vn, inn, P, dvs, dis] = pas_pinn(obj, v, i, s, vin)
            % Un pas Tc du PINN a charge nominale, avec ses derivees par
            % rapport a (v, iL) et a s en mode direct (memes formules que
            % pinn_pas dans banc_pinn_pid.py).
            z = [1.0, (v - 100.0) / 60.0, (i - 12.0) / 10.0, 2 * s - 1, (vin - 195.0) / 45.0, ...
                 (obj.GNOM - obj.GC) / obj.GE];
            T = zeros(3, 6); T(1, 2) = 1 / 60.0; T(2, 3) = 1 / 10.0; T(3, 4) = 2.0;   % dz / d(v, iL, s)
            h = z;
            for l = 1:obj.nc - 1
                h = tanh(h * obj.Wc{l} + obj.bc{l});
                T = (T * obj.Wc{l}) .* (1 - h .* h);
            end
            N = h * obj.Wc{obj.nc} + obj.bc{obj.nc};
            TN = T * obj.Wc{obj.nc};
            vn = v + obj.SVP * N(1); inn = i + obj.SIP * N(2);
            P = [1 + obj.SVP * TN(1,1), obj.SVP * TN(2,1); obj.SIP * TN(1,2), 1 + obj.SIP * TN(2,2)];
            dvs = obj.SVP * TN(3,1); dis = obj.SIP * TN(3,2);
        end

        function [J, gx] = horizon(obj, x, x0, bloc0, r, k1, vin)
            % Simule NH pas a partir de l'etat x0 = [v iL] (instant k1) avec
            % le bloc PID parallele de gains K_DEPART .* x, partant des etats
            % du bloc bloc0 = [xI xF], et le predicteur. Renvoie le cout
            %   J = moyenne de 1/2 (e^2 + W_U u^2) + W_F |x - CENTRE|^2
            % et son gradient par rapport a x (retropropagation dans le
            % temps). Copie de horizon() du banc.
            NHh = obj.NH;
            WUh = obj.WU;
            pente = obj.PSAT;
            Kh = obj.K0 .* x;
            Pg = Kh(1); Ig = Kh(2); Dg = Kh(3);
            xIh = bloc0(1); xFh = bloc0(2);
            v = x0(1); i = x0(2);
            pe = zeros(NHh, 1); pu = pe; pdu = pe; pds = pe; pdvs = pe; pdis = pe;
            pint = false(NHh, 1); pP = zeros(NHh, 4);
            J = 0;
            for j = 1:NHh
                e = r - v;                            % erreur predite
                derivee = obj.NFILTRE * (Dg * e - xFh);   % bloc PID : memes operations que BlocPID
                b = Pg * e + xIh + derivee;
                if b < obj.DMIN
                    u = obj.DMIN; du = pente;
                elseif b > obj.DMAX
                    u = obj.DMAX; du = pente;
                else
                    u = b; du = 1.0;
                end
                entree_I = Ig * e;
                integre = ~((b ~= u) && (sign(entree_I) == sign(b - u)));   % clamping
                ph = mod((k1 + j - 1) * obj.NREG, obj.NPER);   % phase de la porteuse
                y = u * obj.NPER - ph;                % fraction de conduction continue et sa derivee
                if y <= 0
                    s = 0.0; ds = 0.0;
                elseif y >= obj.NREG
                    s = 1.0; ds = 0.0;
                else
                    s = y / obj.NREG; ds = obj.NPER / obj.NREG;
                end
                if obj.PRED_PINN
                    [vn, inn, Pm, dvs, dis] = pas_pinn(obj, v, i, s, vin);
                else
                    [vn, inn, Pm, dvs, dis] = pas_modele(obj, v, i, s, vin);
                end
                J = J + 0.5 * (e * e + WUh * u * u);
                pe(j) = e; pu(j) = u; pdu(j) = du; pds(j) = ds; pint(j) = integre;
                pP(j, :) = [Pm(1,1), Pm(1,2), Pm(2,1), Pm(2,2)];
                pdvs(j) = dvs; pdis(j) = dis;
                if integre
                    xIh = xIh + obj.TC * entree_I;
                end
                xFh = xFh + obj.TC * derivee;
                v = vn; i = inn;
            end
            dxc = x - obj.CENTRE;                     % ecart au centre de la regularisation
            J = J / NHh + obj.WF * (dxc * dxc');
            gK = [0, 0, 0];                           % gradient par rapport a [P I D]
            lv = 0; li = 0; lxI = 0; lxF = 0;         % adjoints de v, iL, xI, xF (etat suivant)
            for j = NHh:-1:1
                e = pe(j); u = pu(j);
                gu = WUh * u / NHh + (lv * pdvs(j) + li * pdis(j)) * pds(j);
                gb = gu * pdu(j);
                gder = gb + lxF * obj.TC;             % derivee entre dans u et dans xF'
                gK(1) = gK(1) + gb * e;
                gK(3) = gK(3) + gder * obj.NFILTRE * e;
                ge = e / NHh + gb * Pg + gder * obj.NFILTRE * Dg;
                if pint(j)
                    gK(2) = gK(2) + lxI * obj.TC * e;
                    ge = ge + lxI * obj.TC * Ig;
                end
                gxI = gb + lxI;
                gxF = lxF - gder * obj.NFILTRE;
                lv_new = pP(j, 1) * lv + pP(j, 3) * li - ge;
                li_new = pP(j, 2) * lv + pP(j, 4) * li;
                lv = lv_new; li = li_new; lxI = gxI; lxF = gxF;
            end
            gx = gK .* obj.K0 + 2 * obj.WF * dxc;
        end

        function J = optimiser(obj, r)
            % N_ADAM iterations d'Adam sur x depuis x courant (depart a chaud),
            % horizon depuis l'etat estime et les etats du bloc PID,
            % projection sur la boite apres chaque iteration.
            x0 = obj.xo(1:2); vin = obj.xo(3);
            bloc0 = [obj.xI, obj.xF];
            th = obj.x;
            m = [0, 0, 0]; vv = [0, 0, 0];            % moments d'Adam, remis a zero (Ito et Wasa, p. 5)
            b1t = 1.0; b2t = 1.0;                     % BETA1^it et BETA2^it, par produits
            J = NaN;
            for it = 1:obj.NADAM
                b1t = b1t * obj.BETA1; b2t = b2t * obj.BETA2;
                [J, g] = horizon(obj, th, x0, bloc0, r, obj.k, vin);
                m = obj.BETA1 * m + (1 - obj.BETA1) * g;
                vv = obj.BETA2 * vv + (1 - obj.BETA2) * g .* g;
                pas_adam = -obj.ALPHA * (m / (1 - b1t)) ./ (sqrt(vv / (1 - b2t)) + obj.EPSA);
                th = min(max(th + pas_adam, obj.XMIN), obj.XMAX);   % projection sur la boite
            end
            obj.x = th;
            obj.K = obj.K0 .* th;
        end
    end

    methods (Static, Access = private)
        function s = charger(nom)
            % Charge un fichier .mat du dossier courant, avec un message clair s'il manque.
            if exist(nom, 'file') ~= 2
                error('pinn_pid_adaptatif : %s introuvable (il doit etre dans le dossier du modele).', nom);
            end
            s = load(nom);
        end
    end
end
