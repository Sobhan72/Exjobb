classdef Solver
    properties
        edof; ex; ey; ed; a; epm
        ndof; nel; bcS; disp
        De; Ds; Dt
        H; sig_y0; sig_eff; r2tol
        P; C
        r1; r1tol; N
        eps; sig; ep
        
        % sig_old = zeros(4,1);
    end

    methods
        function obj = Solver(p)  % Constructor
            % fields = fieldnames(p);
            % for i = 1:numel(fields)
            %     obj.(fields{i}) = p.(fields{i});
            % end

            % Mesh
            [~, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le);
            obj.ndof = 2*((p.lx/p.le + 1)*(p.ly/p.le + 1));
            obj.nel = round(p.lx*p.ly/p.le^2);

            obj.epm = p.epm;
            G12 = p.E/(2*(1+p.v));
            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            obj.Ds = obj.De;
            obj.Dt = repmat(obj.De, 4, 1);
            obj.P = [p.Fco+p.Gco -p.Fco -p.Gco 0; -p.Fco p.Fco+p.Hco -p.Hco 0; -p.Gco -p.Hco p.Gco+p.Hco 0 ; 0 0 0 2*p.Lco];
            obj.H = p.H;
            obj.sig_y0 = p.sig_y0;

            obj.r2tol = p.r2tol;
            obj.r1tol = p.r1tol;
            obj.N = p.N;   
            obj.r1 = ones(obj.ndof, 1);
            obj.disp = p.disp;
            obj.disp(:, 2) = obj.disp(:, 2)/obj.N;
            obj.a = zeros(obj.ndof, 1);
            obj.bcS = obj.addBC(bc, p.ly, p.le, obj.ndof);

            obj.eps = zeros(obj.nel*4, 4);
            obj.sig = zeros(obj.nel*4, 4);
            obj.ep = 0;
            obj.sig_eff = 0;


        end

        function obj = newt(obj)
            for n = 1:obj.N
                bcD = obj.disp;
                Nr = 0;
                while norm(obj.r1) > obj.r1tol || Nr == 0
                    Nr = Nr + 1;
                    obj = FEM(obj, bcD);
                    bcD(:, 2) = bcD(:, 2)*0;
                    fprintf("Nr: %i, r1: %4.2g \n", [Nr, norm(obj.r1)]);
                end
            end
        end

        function obj = FEM(obj, bcD)
            K = zeros(obj.ndof);
            for el = 1:obj.nel
                Ke = plani4e(obj.ex(el,:), obj.ey(el,:), obj.epm, obj.Dt);
                indx = obj.edof(el, 2:end);
                K(indx, indx) = K(indx, indx) + Ke;
            end

            bc = [obj.bcS; bcD];
            da = solveq(K, -obj.r1, bc);
            obj.a = obj.a + da;
            obj.ed = extract_ed(obj.edof, obj.a);

            f_int = zeros(obj.ndof, 1);
            for el = 1:obj.nel
                % fprintf("El: %i \n", el)
                [~, eps2] = plani4s(obj.ex(el,:), obj.ey(el,:), obj.epm, eye(4), obj.ed(el,:));
                for gp = 1:4
                    indxgp = 4*(el-1) + gp;
                    deps = (eps2(gp, :) - obj.eps(indxgp, :))';
                    [obj, obj.sig(indxgp, :), obj.Dt(gp*4-3:gp*4, :)] = hill(obj, deps, eps2(gp, :)', obj.sig(indxgp, :)');
                end
                obj.eps(indxgp-3:indxgp, :) = eps2;
                fe_int = plani4f(obj.ex(el, :), obj.ey(el, :), obj.epm, obj.sig(indxgp-3:indxgp, :))';
                indx = obj.edof(el, 2:end);
                f_int(indx) = f_int(indx) + fe_int;
            end
            f_int(bc(:, 1)) = f_int(bc(:, 1))*0;
            obj.r1 = f_int;
        end

        function [obj, siggp, Dgp] = hill(obj, deps, epsgp, siggp)
            siggp = obj.De*deps + siggp;
            sig_eff_t = sqrt(obj.sig_y0^2*siggp'*obj.P*siggp);
            
            if sig_eff_t > obj.sig_y0
                % e = [epsgp(1:3)-mean(epsgp(1:3)); 2*epsgp(4)];
                [obj, Dgp] = DMat(obj, epsgp);
                siggp = obj.Ds*epsgp;
            else
                obj.sig_eff = sig_eff_t;
                Dgp = obj.De;
            end
            % sig2 = obj.Dgp*deps + obj.sig_old;
            % fprintf("Diff: %5.3g \n", (sig2-siggp)')
            % obj.sig_old = siggp;
        end

        function [obj, Dgp] = DMat(obj, eps)
            epst = eps'*obj.Ds*obj.P*obj.Ds*eps;
            r2 = obj.sig_eff - obj.sig_y0*sqrt(epst);
            iter = 0;
            if isnan(norm(r2))
                error("Residual is NaN")
            end
            while norm(r2) > obj.r2tol
                iter = iter + 1;
                obj.Ds = inv(obj.C + obj.sig_y0^2/obj.sig_eff*obj.ep*obj.P);
                dDsdep = -obj.Ds*obj.P*obj.Ds*(obj.sig_y0^2*(obj.sig_eff-obj.ep*obj.H)/obj.sig_eff^2);
                detdDs = 2*obj.P*obj.Ds*eps*eps';
                epst = eps'*obj.Ds*obj.P*obj.Ds*eps;
                drdep = obj.H - obj.sig_y0/(2*sqrt(epst))*trace(detdDs*dDsdep);
                delta_ep = -r2/drdep;
                obj.ep = obj.ep + delta_ep;
                obj.sig_eff = obj.sig_y0 + obj.H*obj.ep;
                obj.Ds = inv(obj.C + obj.sig_y0^2/obj.sig_eff*obj.ep*obj.P);
                epst = eps'*obj.Ds*obj.P*obj.Ds*eps;
                r2 = obj.sig_eff - obj.sig_y0*sqrt(epst);
                fprintf("  iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            drdeps = -obj.sig_y0*1/sqrt(epst)*obj.Ds*obj.P*obj.Ds*eps;
            depdeps = -drdeps/drdep;
            Dgp = obj.Ds + dDsdep*eps*depdeps';
        end

        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end




% function sig = update_stress(obj, eps)
        %     e = [eps(1:3)-mean(eps(1:3)); 2*eps(4)];
        %     sigkk = 3*obj.Kmod*sum(eps(1:3));
        %     s = obj.G*e;
        %     sig = [s(1:3)+sigkk/3; s(4)];
        % end

        % function sig = update_stress_el(obj, eps)
        %     e = [eps(1:3)-mean(eps(1:3)); eps(4)];
        %     sigkk = 3*obj.Kmod*sum(eps(1:3));
        %     s = 2*obj.Ge*e;
        %     sig = [s(1:3)+sigkk/3; s(4)];
        % end

        % function stress_eff_out = stress_eff(~, sig)    % Calculate effective strain
        %     s = [sig(1:3)-mean(sig(1:3)); sig(4)];
        %     J2 = 0.5*(s'*s + s(4)*s(4));
        %     stress_eff_out = sqrt(J2*3);
        % end
             
        % function D_ep = Dtan(obj, sig, sig_eff)   %%FEL%%
        %     s = [sig(1:3)-mean(sig(1:3)); sig(4)];
        %     A = obj.H + (obj.sig_y0^2/sig_eff)^2*s'*obj.P*D*obj.P*s; %sig_y = sig_eff
        %     D_p = 1/A*(obj.sig_y0^2/sig_eff)^2*D*obj.P*s*s'*obj.P*D; 
        %     D_ep = D-D_p;
        % end