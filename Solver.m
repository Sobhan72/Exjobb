classdef Solver
    properties
        edof; ex; ey
        ndof; nel; bc
        De; Ge; Kmod; ep
        P
        res
    end
    
    methods
        function obj = Solver(p)  % Constructor
            obj.ep = p.ep;
            obj.De = hooke(obj.ep(1), p.E, p.v);
            obj.Kmod = p.E/(3*(1-2*p.v)); 
            obj.Ge = p.E/(2*(1+p.v)); 
            obj.P = [p.Fco+p.Gco -p.Fco -p.Gco 0; -p.Fco p.Fco+p.Hco -p.Hco 0; -p.Gco -p.Hco p.Gco+p.Hco 0 ; 0 0 0 2*p.Lco];


            % Mesh
            [coord, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le);
            obj.ndof = 2*((p.lx/p.le + 1)*(p.ly/p.le + 1));
            obj.nel = round(p.lx*p.ly/p.le^2);
            obj.bc = obj.addBC(bc, p.ly, p.le, obj.ndof);
        end

        function obj = FEM(obj)
            K = zeros(obj.ndof);
            for el = 1:obj.nel
                Ke = plani4e(obj.ex(el,:), obj.ey(el,:), obj.ep, obj.D); 
                indx = obj.edof(:, 2:end);
                K(indx, indx) = K(indx, indx) + Ke;
            end
            
            a = solveq(K, -obj.res, obj.bc);
            ed = extract_ed(obj.edof, a);

            f_int = zeros(obj.ndof, 1);
            for el = 1:obj.nel
                [~, eps] = plani4s(obj.ex(el,:), obj.ey(el,:), obj.ep, obj.D, ed(el,:));
                [G, sig_eff, e_p_eff] = Gp(sig_eff, eps_eff, e_p_eff, sig_y0, G, Ge, H, rtol, P);
                sig = update_stress(obj, eps, G);
                D = Dtan();

                fe_int = plani4f(obj.ex(el, :), obj.ey(el, :), obj.ep, sig);
                indx = obj.edof(el, 2:end);
                f_int(indx) = f_int(indx) + fe_int;
            end


            obj.res = f_int;
        end

        function [G, sig_eff, e_p_eff] = Gp(sig_eff, eps_eff, e_p_eff, sig_y0, G, Ge, H, rtol, P)
            r = sig_eff - 3*G*eps_eff;
            iter = 0;
            T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2];
            while norm(r) > rtol
                iter = iter + 1;
                Gm = inv(1/Ge*eye(4) + 2*sig_y0^2/sig_eff*e_p_eff*P/T);
                dGdepm = -Gm*(2*sig_y0^2/sig_eff*P/T)*Gm;
                dGdep = dGdepm(1);
                drdep = H - dGdep;
                delta_e_p_eff = -r/drdep;
                e_p_eff = e_p_eff + delta_e_p_eff;
                sig_eff = sig_y0 + H*e_p_eff;
                Gm = inv(1/Ge*eye(4) + 2*sig_y0^2/sig_eff*e_p_eff*P/T);
                G = Gm(1);
                r = sig_eff - 3*G*eps_eff;
                % drde = -3*G;
                % dG = -drde/drdep;
                fprintf("  iter: %i, r: %4.2g \n", [iter, norm(r)])
            end
        end

        function sig = update_stress(obj, eps, G)
            e = [eps(1:3)-mean(eps(1:3)); eps(4)];
            sig_kk = 3*obj.Kmod*sum(eps(1:3));
            s = 2*G*e;
            sig = [s(1:3)+sig_kk/3; s(4)];
        end
            
        function D_ep = Dtan(sig, sig_eff)
            s = [sig(1:3)-mean(sig(1:3)); sig(4)];
            A = obj.H + (obj.sig_y0^2/sig_eff)^2*s'*obj.P*obj.D*obj.P*s; %sig_y = sig_eff
            D_p = 1/A*(obj.sig_y0^2/sig_eff)^2*obj.D*obj.P*s*s'*obj.P*obj.D; 
            D_ep = obj.D-D_p;
        end

        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            disp = 2;
            bc = [bc; [fix zeros(2*nR,1); disp 1]];

        end
    end
end

