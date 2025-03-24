classdef Solver
    properties %#ok<*MINV>
        edof; ex; ey; ed; a
        t; ngp
        ndof; nel; bcS; disp
        De; Ds; Dt; X; iX; Gam
        H; sigy0; 
        sigy; r2tol; DP
        P; C
        r1; r1tol; N
        eps; sig; ep %sige
        epsi; sigi; epi; Dsi; sigyi %sigei
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
            
            obj.t = p.t;
            obj.ngp = p.ngp;
            G12 = p.E/(2*(1+p.v));
            obj.P = [p.Fco+p.Gco -p.Fco -p.Gco 0; -p.Fco p.Fco+p.Hco -p.Hco 0; -p.Gco -p.Hco p.Gco+p.Hco 0 ; 0 0 0 2*p.Lco];
            obj.H = p.H;
            obj.sigy0 = p.sigy0;

            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            [obj.X, obj.iX, obj.Gam] = diagDs(obj);
            obj.Ds = repmat(obj.De, obj.nel*obj.ngp, 1);
            obj.Dsi = repmat(obj.De, obj.nel*obj.ngp, 1);
            obj.Dt = repmat(obj.De, obj.nel*obj.ngp, 1);

            obj.r2tol = p.r2tol;
            obj.DP = p.DP;
            obj.r1tol = p.r1tol;
            obj.N = p.N;   
            obj.r1 = sparse(obj.ndof, 1);
            obj.disp = p.disp;
            obj.disp(:, 2) = obj.disp(:, 2)/obj.N;
            obj.a = zeros(obj.ndof, 1);
            obj.bcS = obj.addBC(bc, p.ly, p.le, obj.ndof);

            obj.eps = zeros(obj.nel*obj.ngp, 4);
            obj.sig = zeros(obj.nel*obj.ngp, 4);
            obj.ep = zeros(obj.nel*obj.ngp, 1);
            % obj.sige = zeros(obj.nel*obj.ngp, 1);
            obj.sigy = p.sigy0*ones(obj.nel*obj.ngp, 1);

            obj.epsi = zeros(obj.nel*obj.ngp, 4);
            obj.sigi = zeros(obj.nel*obj.ngp, 4);
            obj.epi = zeros(obj.nel*obj.ngp, 1);
            % obj.sigei = zeros(obj.nel*obj.ngp, 1);
            obj.sigyi = p.sigy0*ones(obj.nel*obj.ngp, 1);
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
                obj.eps = obj.epsi;
                obj.sig = obj.sigi;
                % obj.sige = obj.sigei;
                obj.ep = obj.epi;
                obj.Ds = obj.Dsi;
                obj.sigy = obj.sigyi;
            end
        end

        function obj = FEM(obj, bcD)
            tripK = zeros(obj.nel*64,3);
            for el = 1:obj.nel
                ke = zeros(8);
                for gp = 1:obj.ngp
                    [B, J] = NablaB(obj, gp, el);
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    ke = ke + B'*obj.Dt(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                end
                [rows, cols] = ndgrid(obj.edof(el, :));
                % ke = plani4e(obj.ex(el,:), obj.ey(el,:), obj.epm, obj.Dt(4*obj.ngp*(el-1)+1:4*obj.ngp*el,:));
                tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
            end
            K = sparse(tripK(:,1), tripK(:,2), tripK(:,3), obj.ndof, obj.ndof);

            bc = [obj.bcS; bcD];
            da = obj.solvelin(K, -obj.r1, bc);
            obj.a = obj.a + da;
            obj.ed = obj.a(obj.edof);

            tripf = zeros(obj.nel*8,2);
            for el = 1:obj.nel
                % fprintf("El: %i \n", el)
                % [~, obj.epsi((el-1)*obj.ngp+1:obj.ngp*el,:)] = plani4s(obj.ex(el,:), obj.ey(el,:), obj.epm, eye(4), obj.ed(el,:));
                fein = zeros(8,1);
                for gp = 1:obj.ngp
                    ix = obj.ngp*(el-1) + gp;
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    [B, J] = NablaB(obj, gp, el);
                    obj.epsi(ix, [1 2 4]) = B*obj.ed(el, :)';
                    deps = (obj.epsi(ix, :) - obj.eps(ix, :))';
                    sigtr = obj.De*deps + obj.sig(ix, :)';
                    if sqrt(obj.sigy0^2*sigtr'*obj.P*sigtr) > obj.sigy(ix)
                        % [obj.Dt(ixM, :), obj.sigi(ix, :), obj.Dsi(ixM, :), obj.epi(ix), obj.sigyi(ix)]...
                        % = hill(obj, Deps, obj.epsi(ix, :)', obj.sig(ix, :)', obj.Ds(ixM, :), obj.ep(ix), obj.sigy(ix));
                        if obj.DP
                            [obj.Dt(ixM, :), obj.sigi(ix, :), obj.Dsi(ixM, :), obj.epi(ix)] = DPMat(obj, obj.epsi(ix, :)',obj.Ds(ixM, :), obj.ep(ix));
                        else
                            [obj.Dt(ixM, :), obj.sigi(ix, :), obj.sigyi(ix)] = EPMat(obj, sigtr, obj.sigy(ix));
                        end
                    else
                        obj.sigi(ix, :) = sigtr;
                        % obj.sigei(ix) = sigetr;
                        obj.Dt(ixM, :) = obj.De;
                    end
                    fein = fein + B'*obj.sigi(ix, [1 2 4])'*J*obj.t;
                end
                % fein = plani4f(obj.ex(el, :), obj.ey(el, :), obj.epm, obj.sigi((el-1)*obj.ngp+1:obj.ngp*el, :))';
                tripf((el-1)*8+1:el*8,:) = [obj.edof(el, :)', fein];
            end
            fin = sparse(tripf(:,1), 1, tripf(:,2), obj.ndof, 1);
            fin(bc(:, 1)) = 0;
            obj.r1 = fin;
        end

        % function vararg = hill(obj, deps, epsgp, siggp, Dsgp, epgp, sigygp) %[Dtgp, siggp, Dsgp, epgp, sigygp]
        %     sigtr = obj.De*deps + siggp;
        %     sigetr = sqrt(obj.sigy0^2*sigtr'*obj.P*sigtr);
        %     if sigetr > sigygp
        %         if obj.DP
        %             [vararg{1}, vararg{2}, vararg{3}, vararg{4}] = DPMat(obj, epsgp, Dsgp, epgp);
        %         else
        %             [vararg{1}, vararg{2}, vararg{3}] = EPMat(obj, sigtr, sigygp);
        %         end
        %     else
        %         % sigegp = sigetr;
        %         vararg{1} = obj.De;
        %         vararg{2} = sigtr;
        %         vararg{3} = Dsgp;
        %     end
        % end

        function [Dt, sig, sigy] = EPMat(obj, sigtr, sigy)
            Dep = 0;
            sigt = sigtr'*obj.P*sigtr;
            r2 = sigy - obj.sigy0*sqrt(sigt);
            iter = 0;
            I = eye(4);
            while norm(r2) > obj.r2tol || iter == 0
                iter = iter + 1;
                V = inv(I + obj.sigy0^2*(Dep/(sigy+obj.H*Dep))*obj.Gam);
                U = obj.X*V*obj.iX;
                dUdDep = -obj.X*V*obj.Gam*V*obj.iX*(obj.sigy0^2*sigy/(sigy+obj.H*Dep)^2);
                dsigtdU = 2*obj.P*U*sigtr*sigtr';
                sigt = sigtr'*U*obj.P*U'*sigtr;
                drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dUdDep*dsigtdU);
                DDep = -r2/drdDep;
                Dep = Dep + DDep;
                U = obj.X*inv(I + obj.sigy0^2*(Dep/(sigy+obj.H*Dep))*obj.Gam)*obj.iX;
                sigt = sigtr'*U*obj.P*U'*sigtr;
                r2 = sigy + obj.H*Dep - obj.sigy0 *sqrt(sigt);           
                % fprintf("    iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            sig =  obj.X*inv(I + obj.sigy0^2*(Dep/(sigy+obj.H*Dep))*obj.Gam)*obj.iX*sigtr;
            sigy = sigy + obj.H*Dep;
            dUdDep = -obj.X*V*obj.Gam*V*obj.iX*(obj.sigy0^2*sigy/(sigy+obj.H*Dep)^2);
            dsigtdU = 2*obj.P*U*sigtr*sigtr';
            drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dUdDep*dsigtdU);
            drdDeps = -obj.sigy0/sqrt(sigt)*obj.De*U*obj.P*U*sigtr;
            dDepdDeps = -drdDeps/drdDep;
            Dt = U*obj.De + dUdDep*sigtr*dDepdDeps';
        end

        function [Dt, sig, Ds, ep] = DPMat(obj, eps, Ds, ep)
            epst = eps'*Ds*obj.P*Ds*eps;
            sige = obj.sigy0 + obj.H*ep;
            r2 = sige - obj.sigy0*sqrt(epst);
            iter = 0;
            if isnan(norm(r2))
                error("Residual is NaN")
            end
            while norm(r2) > obj.r2tol || iter == 0
                iter = iter + 1;
                % Ds = inv(obj.C + obj.sigy0^2/sige*ep*obj.P);
                Ds = obj.X*diag(1./diag(eye(4)+obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                dDsdep = -Ds*obj.P*Ds*(obj.sigy0^2*(sige-ep*obj.H)/sige^2);
                detdDs = 2*obj.P*Ds*eps*eps'; 
                epst = eps'*Ds*obj.P*Ds*eps;
                drdep = obj.H - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep);
                Dep = -r2/drdep;
                ep = ep + Dep;
                sige = obj.sigy0 + obj.H*ep;
                % Ds = inv(obj.C + obj.sigy0^2/sige*ep*obj.P);
                Ds = obj.X*diag(1./diag(eye(4)+obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                epst = eps'*Ds*obj.P*Ds*eps;
                r2 = sige - obj.sigy0*sqrt(epst);
                % fprintf("    iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            sig = Ds*eps;
            detdDs = 2*obj.P*Ds*eps*eps'; 
            dDsdep = -Ds*obj.P*Ds*(obj.sigy0^2*(sige-ep*obj.H)/sige^2);
            drdep = obj.H - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep);
            drdeps = -obj.sigy0/sqrt(epst)*Ds*obj.P*Ds*eps;
            depdeps = -drdeps/drdep;
            Dt = Ds + dDsdep*eps*depdeps';
        end

        function [Bgp, detJ] = NablaB(obj, gp, el) % Spatial gradient of shape functions
            node = sqrt(1/3);
            gps = [-node -node; -node node; node node; node -node];
            eta = gps(gp,1); xi = gps(gp,2);
            Npz = 1/4*[eta-1, 1-eta, 1+eta, -1-eta;
                       xi-1, -1-xi, 1+xi, 1-xi];
            Jt = [Npz*obj.ex(el,:)', Npz*obj.ey(el,:)'];
            detJ = det(Jt);
            Npx = Jt\Npz;
            Bgp = zeros(3,8);
            Bgp(1,1:2:8) = Npx(1,:);
            Bgp(2,2:2:8) = Npx(2,:);
            Npx_f = flip(Npx);
            Bgp(3,:) = Npx_f(:);
        end

        function a = solvelin(~,K,f,bc)
            % PURPOSE
            %  Solve static FE-equations considering boundary conditions
            %  for multiple load vectors.
            %
            % INPUT: K : global stiffness matrix, dim(K) = nd x nd
            %        f : global load vector, dim(f) = nd x n      
            %        bc : boundary condition matrix
            %            dim(bc) = nbc x 2 
            %
            % OUTPUT:  a : solution including boundary values
            %              dim(a) = nd x n, 
            %          nd : number of dof's
            %          nbc : number of b.c.'s
            %          n : number of load vectors  
            [nr, nc] = size(f);
            if nargin == 2
                a = K\f;
            else
                fdof = (1:nr)';
                fdof(bc(:,1)) = [];
                s = K(fdof,fdof)\(f(fdof,:)-K(fdof,bc(:,1))*bc(:,2));
                a = zeros([nr, nc]);
                a(fdof,:) = s;
                a(bc(:,1),:) = repmat(bc(:,2),1,nc);
            end
        end

        function [X, iX, Gam] = diagDs(obj)
            [Q, L] = eig(obj.De);
            sL = sqrt(L);
            B = sL*Q'*obj.P*Q*sL;
            [R, Gam] = svd(B);
            X = Q*sL*R;
            iX = inv(X);
        end
        
        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end