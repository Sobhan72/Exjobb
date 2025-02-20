classdef Solver
    properties
        edof; ex; ey; ed
        ndof; nel; bcS
        De; Dp; Kmod; epm
        H; sig_y0; sig_eff
        P; T; C
        res; rtol
        eps; sig; ep
    end

    methods
        function obj = Solver(p)  % Constructor
            obj.epm = p.epm;
            G12 = p.E/(2*(1+p.v)); 
            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            obj.Dp = obj.De;
            obj.Kmod = p.E/(3*(1-2*p.v)); 
            obj.P = [p.Fco+p.Gco -p.Fco -p.Gco 0; -p.Fco p.Fco+p.Hco -p.Hco 0; -p.Gco -p.Hco p.Gco+p.Hco 0 ; 0 0 0 2*p.Lco];
            obj.H = p.H;
            obj.sig_y0 = p.sig_y0;
            obj.sig_eff = 0;
            obj.rtol = p.rtol;
            obj.T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2];
            % obj.G = inv(1/2/obj.Ge*obj.T);

            % Mesh
            [~, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le);
            obj.ndof = 2*((p.lx/p.le + 1)*(p.ly/p.le + 1));
            obj.nel = round(p.lx*p.ly/p.le^2);
            obj.bcS = obj.addBC(bc, p.ly, p.le, obj.ndof);

            obj.eps = zeros(obj.nel*4, 4);
            obj.sig = zeros(4,4);
            obj.ep = 0;
        end

        function obj = FEM(obj, Dt, bcD)
            K = zeros(obj.ndof);
            bc = [obj.bcS; bcD];
            for el = 1:obj.nel
                Ke = plani4e(obj.ex(el,:), obj.ey(el,:), obj.epm, Dt);
                indx = obj.edof(el, 2:end);
                K(indx, indx) = K(indx, indx) + Ke;
            end
            
            a = solveq(K, -obj.res, bc);
            obj.ed = extract_ed(obj.edof, a);

            f_int = zeros(obj.ndof, 1);
            for el = 1:obj.nel
                [~, eps2] = plani4s(obj.ex(el,:), obj.ey(el,:), obj.epm, Dt, obj.ed(el,:));
                for gp = 1:4
                    indx = 4*(el-1) + gp;
                    deps = (eps2(gp, :) - obj.eps(indx, :))';
                    [obj, obj.sig(gp, :)] = hill(obj, deps, eps2(gp, :)', obj.sig(gp, :)');
                end
                obj.eps(indx-3:indx, :) = eps2;
                fe_int = plani4f(obj.ex(el, :), obj.ey(el, :), obj.epm, obj.sig)';
                indx = obj.edof(el, 2:end);
                f_int(indx) = f_int(indx) + fe_int;
            end
            f_int(bc(:,1)) = 0;
            obj.res = f_int;
        end

        function [obj, siggp] = hill(obj, deps, epsgp, siggp)
            siggp = obj.De*deps + siggp;
            sig_eff_t = obj.stress_eff(siggp);
            
            if sig_eff_t > obj.sig_y0
                % e = [epsgp(1:3)-mean(epsgp(1:3)); 2*epsgp(4)];
                obj = Ds(obj, epsgp);
                siggp = obj.Dp*epsgp;
            else
                obj.sig_eff = sig_eff_t;
            end
        end

        function obj = Ds(obj, eps)
            epst = eps'*obj.Dp*obj.P*obj.Dp*eps;
            r = obj.sig_eff - obj.sig_y0*sqrt(epst);
            iter = 0;
            while norm(r) > obj.rtol
                iter = iter + 1;
                obj.Dp = inv(obj.C + obj.sig_y0^2/obj.sig_eff*obj.ep*obj.P);
                dDpdep = -obj.Dp*obj.P*obj.Dp*(obj.sig_y0^2*(obj.sig_eff-obj.ep*obj.H)/obj.sig_eff^2);
                detdD = 2*obj.P*obj.Dp*eps*eps';
                epst = eps'*obj.Dp*obj.P*obj.Dp*eps;
                drdep = obj.H - obj.sig_y0/(2*sqrt(epst))*trace(detdD*dDpdep);
                delta_ep = -r/drdep;
                obj.ep = obj.ep + delta_ep;
                obj.sig_eff = obj.sig_y0 + obj.H*obj.ep;
                obj.Dp = inv(obj.C + obj.sig_y0^2/obj.sig_eff*obj.ep*obj.P);
                epst = eps'*obj.Dp*obj.P*obj.Dp*eps;
                r = obj.sig_eff - obj.sig_y0*sqrt(epst);
                fprintf("  iter: %i, r: %4.2g \n", [iter, norm(r)])
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

        function stress_eff_out = stress_eff(~, sig)    % Calculate effective strain
            s = [sig(1:3)-mean(sig(1:3)); sig(4)];
            J2 = 0.5*(s'*s + s(4)*s(4));
            stress_eff_out = sqrt(J2*3);
        end
             
        % function D_ep = Dtan(obj, sig, sig_eff)   %%FEL%%
        %     s = [sig(1:3)-mean(sig(1:3)); sig(4)];
        %     A = obj.H + (obj.sig_y0^2/sig_eff)^2*s'*obj.P*D*obj.P*s; %sig_y = sig_eff
        %     D_p = 1/A*(obj.sig_y0^2/sig_eff)^2*D*obj.P*s*s'*obj.P*D; 
        %     D_ep = D-D_p;
        % end

        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end