classdef Solver
    properties %#ok<*MINV>
        edof; ex; ey; ed; a
        epm; ngp
        ndof; nel; bcS; disp
        De; Ds; Dt
        H; sig_y0; r2tol
        P; C
        r1; r1tol; N
        eps; sig; sige; ep    
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

            obj.epm = [p.ptype p.t p.ir];
            obj.ngp = p.ir^2;
            G12 = p.E/(2*(1+p.v));
            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            obj.Ds = repmat(obj.De, obj.nel*obj.ngp, 1);
            obj.Dt = repmat(obj.De, obj.nel*obj.ngp, 1);
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

            obj.eps = zeros(obj.nel*obj.ngp, 4);
            obj.sig = zeros(obj.nel*obj.ngp, 4);
            obj.ep = zeros(obj.nel*obj.ngp, 1);
            obj.sige = zeros(obj.nel*obj.ngp, 1);
        end

        function obj = newt(obj)
            for n = 1:obj.N
                fprintf("Load step: %i \n", n);
                bcD = obj.disp;
                Nr = 0;
                while norm(obj.r1) > obj.r1tol || Nr == 0
                    Nr = Nr + 1;
                    obj = FEM(obj, bcD);
                    bcD(:, 2) = bcD(:, 2)*0;
                    fprintf("  Nr: %i, r1: %4.2g \n", [Nr, norm(obj.r1)]);
                end
            end
        end

        function obj = FEM(obj, bcD)
            K = zeros(obj.ndof);
            for el = 1:obj.nel
                Ke = plani4e(obj.ex(el,:), obj.ey(el,:), obj.epm, obj.Dt(4*obj.ngp*(el-1)+1:4*obj.ngp*el,:));
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
                for gp = 1:obj.ngp
                    indxgp = obj.ngp*(el-1) + gp;
                    indxMgp = 4*obj.ngp*(el-1)+(gp-1)*4+1:4*obj.ngp*(el-1)+gp*4;
                    deps = (eps2(gp, :) - obj.eps(indxgp, :))';
                    [obj.sig(indxgp, :), obj.Dt(indxMgp, :), obj.sige(indxgp), obj.Ds(indxMgp, :), obj.ep(indxgp)]...
                    = hill(obj, deps, eps2(gp, :)', obj.sig(indxgp, :)', obj.sige(indxgp), obj.Ds(indxMgp, :), obj.ep(indxgp));
                end
                obj.eps(indxgp-obj.ngp+1:indxgp, :) = eps2;
                fe_int = plani4f(obj.ex(el, :), obj.ey(el, :), obj.epm, obj.sig(indxgp-obj.ngp+1:indxgp, :))';
                indx = obj.edof(el, 2:end);
                f_int(indx) = f_int(indx) + fe_int;
            end
            f_int(bc(:, 1)) = f_int(bc(:, 1))*0;
            obj.r1 = f_int;
        end

        function [siggp, Dtgp, sigegp, Dsgp, epgp, Dt2] = hill(obj, deps, epsgp, siggp, sigegp, Dsgp, epgp)
            siggp = obj.De*deps + siggp;
            sig_eff_t = sqrt(obj.sig_y0^2*siggp'*obj.P*siggp);
            
            if sig_eff_t > obj.sig_y0
                [Dt2, sigegp, Dsgp, epgp, Dtgp] = DMat(obj, epsgp, sigegp, Dsgp, epgp);
                siggp = Dsgp*epsgp;
            else
                sigegp = sig_eff_t;
                Dtgp = obj.De;
            end
        end

        function [Dt, sige, Ds, ep, Dt2] = DMat(obj, eps, sige, Ds, ep)
            epst = eps'*Ds*obj.P*Ds*eps;
            r2 = sige - obj.sig_y0*sqrt(epst);
            iter = 0;
            if isnan(norm(r2))
                error("Residual is NaN")
            end
            while norm(r2) > obj.r2tol || iter == 0
                iter = iter + 1;
                Ds = inv(obj.C + obj.sig_y0^2/sige*ep*obj.P);
                dDsdep = -Ds*obj.P*Ds*(obj.sig_y0^2*(sige-ep*obj.H)/sige^2);
                detdDs = 2*obj.P*Ds*eps*eps'; 
                epst = eps'*Ds*obj.P*Ds*eps;
                drdep = obj.H - obj.sig_y0/(2*sqrt(epst))*trace(detdDs*dDsdep);
                delta_ep = -r2/drdep;
                ep = ep + delta_ep;
                sige = obj.sig_y0 + obj.H*ep;
                Ds = inv(obj.C + obj.sig_y0^2/sige*ep*obj.P);
                epst = eps'*Ds*obj.P*Ds*eps;
                r2 = sige - obj.sig_y0*sqrt(epst);
                fprintf("    iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            drdeps = -obj.sig_y0/sqrt(epst)*Ds*obj.P*Ds*eps;
            depdeps = -drdeps/drdep;
            Dt = Ds + dDsdep*eps*depdeps';
            % Dt2 = newDt(obj, sige, Ds, eps, delta_ep);
            Dt2 = Dloop(obj, Ds, dDsdep, depdeps, eps);
        end

        function Dt = newDt(obj, sige, Ds, eps, delta_ep)
            sig = Ds*eps;
            dfdsig = obj.sig_y0^2/sige*obj.P*sig;
            df2dsig2 = obj.sig_y0^2/sige*obj.P*(eye(4)-obj.sig_y0^2/sige^2*sig*sig'*obj.P);
            Da = inv(obj.C + delta_ep*df2dsig2);
            A = dfdsig'*Da*dfdsig + obj.H;
            Dt = Da - 1/A*Da*(dfdsig*dfdsig')*Da;
        end

        function D = Dloop(obj, Ds, dDsdep, depdeps, eps)
            indx = [1 1 1 1; 
                    1 1 2 2; 
                    1 1 3 3; 
                    1 1 1 2; 
                    2 2 2 2; 
                    2 2 3 3;
                    2 2 1 2;
                    3 3 3 3;
                    3 3 1 2;
                    1 2 1 2];

            Dt = [];

            eps = [eps(1) eps(4)/2 0;
                   eps(4)/2 eps(2) 0;
                   0 0 eps(3)];

            depdeps = [depdeps(1) depdeps(4)*2 0;
                       depdeps(4)*2 depdeps(2) 0;
                       0 0 depdeps(3)];
          

            dDsdep = tensor(obj,dDsdep);

            for ii = indx'
                Dti = 0;
                for k = 1:3
                    for l = 1:3
                        Dti = Dti + dDsdep(ii(1), ii(2),k, l)*depdeps(ii(3),ii(4))*eps(k,l);
                    end
                end
                Dt = [Dt, Dti];
            end

            D = Ds + [Dt(1), Dt(2), Dt(3), Dt(4);
                      Dt(2), Dt(5), Dt(6), Dt(7);
                      Dt(3), Dt(6), Dt(8), Dt(9);
                      Dt(4), Dt(7), Dt(9), Dt(10)];
        end

        function t = tensor(obj, matrix)
            t = zeros(3,3,3,3);
            t(1,1,1,1) = matrix(1,1);
            t(1,1,1,2) = matrix(1,4);
            t(1,1,2,1) = matrix(1,4);
            t(1,1,2,2) = matrix(1,2);
            t(1,1,3,3) = matrix(1,3);

            t(1,2,1,1) = matrix(1,4);
            t(1,2,1,2) = matrix(4,4);
            t(1,2,2,1) = matrix(4,4);
            t(1,2,2,2) = matrix(2,4);
            t(1,2,3,3) = matrix(3,4);

            t(2,1,1,1) = matrix(1,4);
            t(2,1,1,2) = matrix(4,4);
            t(2,1,2,1) = matrix(4,4);
            t(2,1,2,2) = matrix(2,4);
            t(2,1,3,3) = matrix(3,4);

            t(2,2,1,1) = matrix(2,1);
            t(2,2,1,2) = matrix(2,4);
            t(2,2,2,1) = matrix(2,4);
            t(2,2,2,2) = matrix(2,2);
            t(2,2,3,3) = matrix(2,3);

            t(3,3,1,1) = matrix(3,1);
            t(3,3,1,2) = matrix(3,4);
            t(3,3,2,1) = matrix(3,4);
            t(3,3,2,2) = matrix(3,2);
            t(3,3,3,3) = matrix(3,3);

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