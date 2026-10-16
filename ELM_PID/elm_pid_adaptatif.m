classdef elm_pid_adaptatif < matlab.System
    %% ELM_PID_ADAPTATIF -- adaptateur ELM des gains du PID
    %
    % CE QUE FAIT CE BLOC
    % -------------------
    % Toutes les 0.5 ms, il ajuste les gains du PID par un gradient
    % normalise. Le sens et l'amplitude de chaque ajustement viennent du
    % jacobien d'un modele ELM du convertisseur, mis a jour en ligne
    % (OS-ELM). C'est la copie exacte de la classe AdaptateurELM de
    % banc_elm_pid.py : memes reglages, memes noms, meme ordre des
    % operations. Tester_ELM_PID_Rejeu.m verifie que les deux donnent les
    % memes gains pas par pas.
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
    % 1. Sensibilite de u(k) a un changement des gains au debut de la
    %    fenetre (derivee exacte de la loi du bloc PID) :
    %      du/dP = e(k) ; du/dI = Tc x (somme des e de la fenetre avant k) ;
    %      du/dD = N (e(k) - phiD), puis phiD = phiD + Tc N (e(k) - phiD).
    % 2. Fenetre de 110 periodes (0.5 ms) : moyennes ybar, dbar, ebar, et
    %    des trois sensibilites ; u en butee (<= D_MIN ou >= D_MAX) rend la
    %    fenetre suspecte.
    % 3. Fin de fenetre : l'ELM predit ybar(n) depuis [ybar(n-1) ybar(n-2)
    %    dbar(n) dbar(n-1) dbar(n-2)] et donne J = d ybar(n) / d dbar(n).
    % 4. Porte : ni la fenetre ni les deux precedentes suspectes.
    %    a. OS-ELM si APPRENDRE et |ybar - prediction| > ZONE_MORTE_ESTIMATION ;
    %    b. si |ebar| > ZONE_MORTE et J > 0 : x = K ./ K_DEPART,
    %       phi = J (s .* K_DEPART) ;
    %       dx = ETA ebar phi / (EPS_PHI + phi phi') + ALPHA (x - x_prec) ;
    %       projection de x + dx sur l'ensemble admissible
    %       (ensemble_gains_elm.mat) ; K = K_DEPART .* x.
    %
    % FICHIERS LUS (dans le dossier courant, un seul exemplaire)
    %   elm_pid_modele.mat      : modele ELM (entrainement_elm.py)
    %   elm_pid_reglages.mat    : reglages (banc_elm_pid.py)
    %   ensemble_gains_elm.mat  : table de l'ensemble admissible
    %                             (ensemble_gains_elm.py)
    % Aucun nombre de l'algorithme n'est ecrit dans ce fichier.
    %
    % GainsFixes = true : le bloc sort les gains de Ziegler-Nichols a chaque
    % pas (pas d'adaptation) ; le modele doit alors redonner le PID
    % classique. Sert seulement aux auto-tests T0 des scripts de
    % construction.
    %
    % Utilise dans Simulink en execution interpretee (Construction_ELM_PID.m).
    % Compatible MATLAB R2024a.

    properties (Nontunable)
        FichierModele = 'elm_pid_modele.mat';
        FichierReglages = 'elm_pid_reglages.mat';
        FichierEnsemble = 'ensemble_gains_elm.mat';
        AfficherBilan (1,1) logical = true;
        GainsFixes (1,1) logical = false;            % true : gains de Ziegler-Nichols (auto-tests T0)
    end

    properties (Access = private)
        W; B; BETA0; P0; XM; XE; TM; TE_; n
        BETA; P
        K0; ETA; ZM; ZME; ALPHA; EPS; LAMBDA; APPRENDRE; JC; GLISSEMENT
        XMIN; XMAX; NBIS; HGRAD; RESTAURATION; NPROJ
        GT; L2MIN; PASL2; NT
        NF; TC; NFILTRE; TCN; DMIN; DMAX
        K; x; x_prec
        cpt; sy; sd; se; hors; SE; phiD; sS
        yb1; yb2; db1; db2; sus1; sus2; nfen
        der_porte; der_adapte; nb_portes; nb_adapt
    end

    methods
        function obj = elm_pid_adaptatif(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function info = lire_fenetre(obj)
            % Etat de la derniere fenetre terminee (pour le test de rejeu).
            info = struct('nfen', obj.nfen, 'porte', obj.der_porte, 'adapte', obj.der_adapte);
        end
    end

    methods (Access = protected)

        function setupImpl(obj)
            m = elm_pid_adaptatif.charger(obj.FichierModele);
            r = elm_pid_adaptatif.charger(obj.FichierReglages);
            t = elm_pid_adaptatif.charger(obj.FichierEnsemble);
            obj.W = double(m.W);
            obj.B = double(m.B(:)');
            obj.BETA0 = double(m.BETA(:));
            obj.P0 = double(m.P0);
            obj.XM = double(m.X_MOY(:)');
            obj.XE = double(m.X_EC(:)');
            obj.TM = double(m.T_MOY);
            obj.TE_ = double(m.T_EC);
            obj.n = size(obj.W, 2);
            if size(obj.W, 1) ~= 5 || numel(obj.BETA0) ~= obj.n + 6 || ~isequal(size(obj.P0), [obj.n + 6, obj.n + 6])
                error('elm_pid_adaptatif : tailles inattendues dans %s.', obj.FichierModele);
            end
            obj.K0 = double(r.K_DEPART(:)');
            obj.ETA = double(r.ETA);
            obj.ZM = double(r.ZONE_MORTE);
            obj.ZME = double(r.ZONE_MORTE_ESTIMATION);
            obj.ALPHA = double(r.ALPHA);
            obj.EPS = double(r.EPS_PHI);
            obj.LAMBDA = double(r.LAMBDA);
            obj.APPRENDRE = double(r.APPRENDRE) ~= 0;
            obj.JC = double(r.JACOBIEN_CONSTANT);
            obj.GLISSEMENT = double(r.GLISSEMENT) ~= 0;
            obj.XMIN = double(r.MULT_MIN);
            obj.XMAX = double(r.MULT_MAX);
            obj.NBIS = double(r.N_BISSECTIONS);
            obj.HGRAD = double(r.H_GRADIENT);
            if ~isfield(r, 'FORMAT') || double(r.FORMAT) ~= 4
                error(['elm_pid_adaptatif : %s n''a pas le format ecrit par banc_elm_pid.py. ' ...
                       'Ne pas melanger les fichiers de ELM_PID_INCREMENTAL et de ELM_PID.'], obj.FichierReglages);
            end
            obj.RESTAURATION = double(r.RESTAURATION) ~= 0;
            obj.NPROJ = double(r.N_PROJECTION);
            obj.NF = double(r.NF);
            obj.TC = double(r.TC);
            obj.NFILTRE = double(r.N_FILTRE);
            obj.TCN = obj.TC * obj.NFILTRE;
            obj.DMIN = double(r.D_MIN);
            obj.DMAX = double(r.D_MAX);
            obj.GT = double(t.G_TABLE);
            obj.L2MIN = double(t.LOG2_MIN);
            obj.PASL2 = double(t.PAS_LOG2);
            obj.NT = size(obj.GT, 1);
            if ~isfield(t, 'FORMAT') || double(t.FORMAT) ~= 2
                error('elm_pid_adaptatif : %s n''a pas le format ecrit par ensemble_gains_elm.py.', ...
                      obj.FichierEnsemble);
            end
            if max(abs(double(t.K_DEPART(:)') - obj.K0) ./ obj.K0) > 1e-12
                error('elm_pid_adaptatif : %s et %s n''ont pas les memes gains de depart.', ...
                      obj.FichierEnsemble, obj.FichierReglages);
            end
            if double(m.NF) ~= obj.NF || abs(double(m.TC) - double(r.TC)) > 1e-18
                error('elm_pid_adaptatif : %s et %s n''ont pas les memes fenetres (NF ou TC).', ...
                      obj.FichierModele, obj.FichierReglages);
            end
        end

        function resetImpl(obj)
            obj.BETA = obj.BETA0;
            obj.P = obj.P0;
            obj.x = [1 1 1];                          % multiplicateurs de Ziegler-Nichols
            obj.x_prec = obj.x;
            obj.K = obj.K0 .* obj.x;
            obj.cpt = 0; obj.sy = 0; obj.sd = 0; obj.se = 0; obj.hors = false;
            obj.SE = 0; obj.phiD = 0; obj.sS = [0 0 0];
            obj.yb1 = 0; obj.yb2 = 0; obj.db1 = 0; obj.db2 = 0;
            obj.sus1 = true; obj.sus2 = true;
            obj.nfen = 0;
            obj.der_porte = false; obj.der_adapte = false; obj.nb_portes = 0; obj.nb_adapt = 0;
        end

        function K = outputImpl(obj, ~, ~, ~)
            % Gains du pas, d'apres l'etat seulement (pas de traversee directe).
            K = obj.K;
        end

        function updateImpl(obj, e, mesure, u)
            % Sensibilites (memes operations, meme ordre que mise_a_jour du banc)
            obj.sS(1) = obj.sS(1) + e;
            obj.sS(2) = obj.sS(2) + obj.TC * obj.SE;
            obj.sS(3) = obj.sS(3) + obj.NFILTRE * (e - obj.phiD);
            obj.SE = obj.SE + e;
            obj.phiD = obj.phiD + obj.TCN * (e - obj.phiD);
            obj.sy = obj.sy + mesure;
            obj.sd = obj.sd + u;
            obj.se = obj.se + e;
            if u <= obj.DMIN || u >= obj.DMAX           % sortie du PID en butee : fenetre suspecte
                obj.hors = true;
            end
            obj.cpt = obj.cpt + 1;
            if obj.cpt >= obj.NF
                fin_de_fenetre(obj);
            end
        end

        function [f1, f2, f3] = isInputDirectFeedthroughImpl(~, ~, ~, ~)
            f1 = false;
            f2 = false;
            f3 = false;
        end

        function releaseImpl(obj)
            if obj.AfficherBilan && ~isempty(obj.nfen) && obj.nfen > 0
                fprintf(['  ELM-PID : %d fenetres, porte ouverte sur %d, pas calcules sur %d ; ' ...
                         'gains finaux P=%.5f I=%.4f D=%.4e (x ZN %.3f %.3f %.3f)\n'], obj.nfen, obj.nb_portes, ...
                        obj.nb_adapt, obj.K(1), obj.K(2), obj.K(3), obj.x(1), obj.x(2), obj.x(3));
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

        function fin_de_fenetre(obj)
            yb = obj.sy / obj.NF;
            db = obj.sd / obj.NF;
            eb = obj.se / obj.NF;
            s = obj.sS / obj.NF;
            sus = obj.hors;
            porte = false;
            adapte = false;
            if obj.nfen >= 2
                xin = [obj.yb1, obj.yb2, db, obj.db1, obj.db2];
                [yhat, J, h] = elm_pid_adaptatif.evaluer_modele(obj.W, obj.B, obj.BETA, ...
                    obj.XM, obj.XE, obj.TM, obj.TE_, xin);
                if ~isnan(obj.JC)
                    J = obj.JC;                       % ablation : jacobien constant
                end
                err = yb - yhat;
                porte = ~sus && ~obj.sus1 && ~obj.sus2;   % porte : saturation seulement
                if porte
                    if obj.APPRENDRE && abs(err) > obj.ZME    % a. OS-ELM, zone morte d'estimation
                        hc = h';
                        Ph = obj.P * hc;
                        gain = Ph / (obj.LAMBDA + hc' * Ph);
                        obj.BETA = obj.BETA + gain * (err / obj.TE_);
                        obj.P = (obj.P - gain * Ph') / obj.LAMBDA;
                    end
                    ancien = obj.x;
                    if ~obj.GainsFixes && abs(eb) > obj.ZM && J > 0   % b. gains
                        phi = J * (s .* obj.K0);
                        dx = obj.ETA * eb * phi / (obj.EPS + phi * phi') + obj.ALPHA * (obj.x - obj.x_prec);
                        obj.x = projeter(obj, obj.x, dx);
                        adapte = true;
                        obj.nb_adapt = obj.nb_adapt + 1;
                    end
                    obj.x_prec = ancien;
                    obj.K = obj.K0 .* obj.x;
                    obj.nb_portes = obj.nb_portes + 1;
                end
            end
            obj.der_porte = porte;
            obj.der_adapte = adapte;
            obj.yb2 = obj.yb1; obj.yb1 = yb;
            obj.db2 = obj.db1; obj.db1 = db;
            obj.sus2 = obj.sus1; obj.sus1 = sus;
            obj.nfen = obj.nfen + 1;
            obj.cpt = 0; obj.sy = 0; obj.sd = 0; obj.se = 0; obj.hors = false;
            obj.SE = 0; obj.phiD = 0; obj.sS = [0 0 0];
        end

        function g = g_interp(obj, x)
            % g de l'ensemble admissible, interpolation trilineaire en log2
            % des multiplicateurs (memes operations que g_interp du banc).
            uu = (log2(x) - obj.L2MIN) / obj.PASL2;
            uu = min(max(uu, 0), obj.NT - 1);
            i = min(floor(uu), obj.NT - 2);
            f = uu - i;
            g = 0;
            for di = 0:1
                if di, wi = f(1); else, wi = 1 - f(1); end
                for dj = 0:1
                    if dj, wj = f(2); else, wj = 1 - f(2); end
                    for dk = 0:1
                        if dk, wk = f(3); else, wk = 1 - f(3); end
                        g = g + wi * wj * wk * obj.GT(i(1) + di + 1, i(2) + dj + 1, i(3) + dk + 1);
                    end
                end
            end
        end

        function gr = gradient_g(obj, x)
            % Gradient de g par differences centrees (normale a la frontiere).
            hh = obj.HGRAD;
            gr = [0 0 0];
            for q = 1:3
                ep = [0 0 0];
                ep(q) = hh;
                gr(q) = (g_interp(obj, x + ep) - g_interp(obj, x - ep)) / (2 * hh);
            end
        end

        function [p, ok] = restaurer(obj, p, gr)
            % Restauration : recul de p le long de -gr jusqu'a la
            % frontiere (doublement puis bissection) ; ok = false si un recul
            % de 1 ne suffit pas (meme algorithme que restaurer() du banc).
            lo_ = obj.XMIN; hi_ = obj.XMAX;
            nh = gr / sqrt(gr * gr');
            s_ok = 1e-6;
            ok = true;
            while g_interp(obj, min(max(p - s_ok * nh, lo_), hi_)) > 0
                s_ok = s_ok * 2;
                if s_ok > 1
                    ok = false;
                    return;
                end
            end
            if s_ok == 1e-6
                s_ko = 0;
            else
                s_ko = 0.5 * s_ok;
            end
            for it = 1:obj.NBIS
                m = 0.5 * (s_ko + s_ok);
                if g_interp(obj, min(max(p - m * nh, lo_), hi_)) <= 0
                    s_ok = m;
                else
                    s_ko = m;
                end
            end
            p = min(max(p - s_ok * nh, lo_), hi_);
        end

        function xn = projeter(obj, xc, dx)
            % Projection du pas sur l'ensemble admissible : point
            % admissible le plus proche de la cible (projection du gradient
            % de Rosen avec restauration) ou, pour l'ablation, glissement
            % puis bissection seuls. Meme algorithme que projeter()
            % du banc.
            lo_ = obj.XMIN; hi_ = obj.XMAX;
            xn = min(max(xc + dx, lo_), hi_);
            if g_interp(obj, xn) <= 0
                return;
            end
            z = xn;
            a = 0; b = 1;
            for it = 1:obj.NBIS
                m = 0.5 * (a + b);
                if g_interp(obj, min(max(xc + m * dx, lo_), hi_)) <= 0
                    a = m;
                else
                    b = m;
                end
            end
            xb = min(max(xc + a * dx, lo_), hi_);
            if ~obj.GLISSEMENT
                xn = xb;
                return;
            end
            if obj.RESTAURATION
                y = xb;
                for it = 1:obj.NPROJ
                    gr = gradient_g(obj, y);
                    if gr * gr' == 0
                        break;
                    end
                    d = z - y;
                    sortant = gr * d';
                    if sortant > 0
                        d = d - sortant / (gr * gr') * gr;
                    end
                    pt = min(max(y + d, lo_), hi_);
                    if g_interp(obj, pt) > 0
                        [pt, ok] = restaurer(obj, pt, gr);
                        if ~ok
                            break;
                        end
                    end
                    if sum((z - pt) .^ 2) >= sum((z - y) .^ 2)
                        break;
                    end
                    y = pt;
                end
                xn = y;
                return;
            end
            gr = gradient_g(obj, xb);
            reste = (1 - a) * dx;
            sortant = gr * reste';
            if sortant > 0
                reste = reste - sortant / (gr * gr') * gr;
            end
            xn = min(max(xb + reste, lo_), hi_);
            if g_interp(obj, xn) <= 0
                return;
            end
            a = 0; b = 1;
            for it = 1:obj.NBIS
                m = 0.5 * (a + b);
                if g_interp(obj, min(max(xb + m * reste, lo_), hi_)) <= 0
                    a = m;
                else
                    b = m;
                end
            end
            xn = min(max(xb + a * reste, lo_), hi_);
        end
    end

    methods (Static)

        function [y, J, h] = evaluer_modele(W, B, BETA, XM, XE, TM, TE_, x)
            % Prediction ybar(n), jacobien d ybar(n) / d dbar(n) et vecteur h
            % (memes formules que AdaptateurELM.modele dans banc_elm_pid.py).
            n = size(W, 2);
            xn = (x - XM) ./ XE;
            g = 1 ./ (1 + exp(-(xn * W + B)));
            h = [g, xn, 1];
            y = (h * BETA) * TE_ + TM;
            dy = ((g .* (1 - g)) .* BETA(1:n)') * W' + BETA(n+1:n+5)';
            J = dy(3) / XE(3) * TE_;
        end
    end

    methods (Static, Access = private)
        function s = charger(nom)
            if exist(nom, 'file') ~= 2
                error('elm_pid_adaptatif : %s introuvable (il doit etre dans le dossier du modele).', nom);
            end
            s = load(nom);
        end
    end
end
