classdef Solver
    properties
        edof; ex; ey; ed; a
        t; ngp
        ndof; nel; bcS; disp
        De; Ds; Dt; X; iX; Gam
        H; sigy0; Kinf; del
        sigy; r2tol; DP
        P; C; DeP
        r1; r1tol; N
        eps; sig; ep
        epsi; sigi; epi; Dsi; sigyi
        Z
    end

    methods
        function obj = Solver(p)  % Constructor
            % Mesh
            [~, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le);
            obj.ndof = 2*((p.lx/p.le + 1)*(p.ly/p.le + 1));
            obj.nel = round(p.lx*p.ly/p.le^2);
            obj.Z = filterMatrix(obj, p.le, p.re);
            
            obj.t = p.t;
            obj.ngp = p.ngp;

            Fco = 1/2*(1/p.sigy01^2+1/p.sigy02^2-1/p.sigy03^2);
            Gco = 1/2*(1/p.sigy01^2-1/p.sigy02^2+1/p.sigy03^2);
            Hco = 1/2*(-1/p.sigy01^2+1/p.sigy02^2+1/p.sigy03^2);
            obj.sigy0 = sqrt(3/(2*(Fco+Gco+Hco)));
            Lco = 3/(2*obj.sigy0^2);
            
            if 4/(p.sigy01^2*p.sigy02^2) <= (1/p.sigy03^2-(1/p.sigy01^2+1/p.sigy02^2))^2
                error("Not positive definite")
            end

            obj.P = [Fco+Gco -Fco -Gco 0; -Fco Fco+Hco -Hco 0; -Gco -Hco Gco+Hco 0; 0 0 0 2*Lco];
            obj.H = p.H;

            obj.Kinf = p.Kinf;
            obj.del = p.del;

            G12 = p.E1/(2*(1+p.v12));
            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            [obj.X, obj.iX, obj.Gam] = diagDs(obj);
            obj.Ds = repmat(obj.De, obj.nel*obj.ngp, 1);
            obj.Dsi = repmat(obj.De, obj.nel*obj.ngp, 1);
            obj.Dt = repmat(obj.De, obj.nel*obj.ngp, 1);

            obj.DeP = obj.De*obj.P;

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
            obj.sigy = p.sigy0*ones(obj.nel*obj.ngp, 1);

            obj.epsi = zeros(obj.nel*obj.ngp, 4);
            obj.sigi = zeros(obj.nel*obj.ngp, 4);
            obj.epi = zeros(obj.nel*obj.ngp, 1);
            obj.sigyi = p.sigy0*ones(obj.nel*obj.ngp, 1);
        end

        %% FEM
        function obj = newt(obj) % Newton-Raphson method
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
                obj.eps = obj.epsi; obj.sig = obj.sigi; obj.ep = obj.epi; obj.Ds = obj.Dsi; obj.sigy = obj.sigyi;
            end
        end

        function obj = FEM(obj, bcD) % Main FEM function
            tripK = zeros(obj.nel*64,3);
            for el = 1:obj.nel
                ke = zeros(8);
                for gp = 1:obj.ngp
                    [B, J] = NablaB(obj, gp, el);
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    ke = ke + B'*obj.Dt(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                end
                [rows, cols] = ndgrid(obj.edof(el, :));
                tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
            end
            K = sparse(tripK(:,1), tripK(:,2), tripK(:,3), obj.ndof, obj.ndof);

            bc = [obj.bcS; bcD];
            da = obj.solvelin(K, -obj.r1, bc);
            obj.a = obj.a + da;
            obj.ed = obj.a(obj.edof);

            tripf = zeros(obj.nel*8,2);
            for el = 1:obj.nel
                fein = zeros(8,1);
                for gp = 1:obj.ngp
                    ix = obj.ngp*(el-1) + gp;
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    [B, J] = NablaB(obj, gp, el);
                    obj.epsi(ix, [1 2 4]) = B*obj.ed(el, :)';
                    deps = (obj.epsi(ix, :) - obj.eps(ix, :))';
                    sigtr = obj.De*deps + obj.sig(ix, :)';
                    if sqrt(obj.sigy0^2*sigtr'*obj.P*sigtr) > obj.sigy(ix)
                        if obj.DP
                            [obj.Dt(ixM, :), obj.sigi(ix, :), obj.Dsi(ixM, :), obj.epi(ix)] = DPMat(obj, obj.epsi(ix, :)',obj.Ds(ixM, :), obj.ep(ix));
                        else
                            [obj.Dt(ixM, :), obj.sigi(ix, :), obj.sigyi(ix), Dep] = EPMat(obj, sigtr, obj.sigy(ix));
                            obj.epi(ix) = obj.epi(ix) + Dep;
                        end
                    else
                        obj.sigi(ix, :) = sigtr;
                        obj.Dt(ixM, :) = obj.De;
                    end
                    fein = fein + B'*obj.sigi(ix, [1 2 4])'*J*obj.t;
                end
                tripf((el-1)*8+1:el*8,:) = [obj.edof(el, :)', fein];
            end
            fin = sparse(tripf(:,1), 1, tripf(:,2), obj.ndof, 1);
            fin(bc(:, 1)) = 0;
            obj.r1 = fin;
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

        function a = solvelin(~,K,f,bc) % Solve FE-equations
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

        %% Material Functions 
        function [Dt, sig, sigy, Dep] = EPMat(obj, sigtr, sigy)
            Dep = 0;
            sigt = sigtr'*obj.P*sigtr;
            r2 = sigy - obj.sigy0*sqrt(sigt);
            iter = 0;
            while norm(r2) > obj.r2tol || iter == 0
                iter = iter + 1;
                V = diag(1./diag(eye(4) + obj.sigy0^2*(Dep/(sigy + obj.H*Dep))*obj.Gam));
                U = obj.X*V*obj.iX;
                %dUdDep = -obj.X*V*obj.Gam*V*obj.iX*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2);
                dUdDep = -U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2);
                dsigtdU = 2*obj.P*U*(sigtr*sigtr');
                sigt = sigtr'*U*obj.P*U'*sigtr;
                drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dUdDep*dsigtdU);
                DDep = -r2/drdDep;
                Dep = Dep + DDep;
                V = diag(1./diag(eye(4) + obj.sigy0^2*(Dep/(sigy + obj.H*Dep))*obj.Gam));
                U = obj.X*V*obj.iX;
                sigt = sigtr'*U*obj.P*U'*sigtr;
                r2 = sigy + obj.H*Dep - obj.sigy0*sqrt(sigt);           
                % fprintf("    iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            sig = U*sigtr;
            sigy = sigy + obj.H*Dep;
            % dUdDep = -obj.X*V*obj.Gam*V*obj.iX*(obj.sigy0^2*sigy/(sigy+obj.H*Dep)^2);
            dUdDep = -U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2);
            dsigtdU = 2*obj.P*U*(sigtr*sigtr');
            drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dUdDep*dsigtdU);
            drdDeps = -obj.sigy0/sqrt(sigt)*obj.De*U'*obj.P*U*sigtr;
            dDepdDeps = -drdDeps/drdDep;
            Dt = U*obj.De + dUdDep*sigtr*dDepdDeps';
        end

        function [Dt, sig, Ds, ep] = DPMat(obj, eps, Ds, ep)
            epst = eps'*Ds*obj.P*Ds*eps;
            sige = obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.del*ep));
            r2 = sige - obj.sigy0*sqrt(epst);
            iter = 0;
            while norm(r2) > obj.r2tol || iter == 0
                iter = iter + 1;
                Ds = obj.X*diag(1./diag(eye(4) + obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                dDsdep = -Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.del*exp(-obj.del*ep))*ep)/sige^2);
                detdDs = 2*obj.P*Ds*(eps*eps');
                epst = eps'*Ds*obj.P*Ds*eps;
                drdep = obj.H + obj.Kinf*obj.del*exp(-obj.del*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep);
                Dep = -r2/drdep;
                ep = ep + Dep;
                sige = obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.del*ep));
                Ds = obj.X*diag(1./diag(eye(4) + obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                epst = eps'*Ds*obj.P*Ds*eps;
                r2 = sige - obj.sigy0*sqrt(epst);
                % fprintf("    iter: %i, r2: %4.2g \n", [iter, norm(r2)])
            end
            sig = Ds*eps;
            detdDs = 2*obj.P*Ds*(eps*eps');
            dDsdep = -Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.del*exp(-obj.del*ep))*ep)/sige^2);
            drdep = obj.H + obj.Kinf*obj.del*exp(-obj.del*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep);
            drdeps = -obj.sigy0/sqrt(epst)*Ds*obj.P*Ds*eps;
            depdeps = -drdeps/drdep;
            Dt = Ds + dDsdep*eps*depdeps';
        end

        function [X, iX, Gam] = diagDs(obj)
            [Q, L] = eig(obj.De);
            sL = sqrt(L);
            B = sL*Q'*obj.P*Q*sL;
            [R, Gam] = svd(B);
            X = Q*sL*R;
            iX = inv(X);
        end

        %% Optimization
        function Z = filterMatrix(obj, le, re)
            ec = [obj.ex(:, 1) + le/2, obj.ey(:, 1) + le/2];
            I = [];
            n = ceil(re/1.4); % 2*n is max kernel size
            [x, y] = meshgrid(-n:n, -n:n);
            weights = max(0, 1 - sqrt(x.^2 + y.^2)/re); % Kernel is hat function
            sw = sum(weights(:));
            r0 = le*re;
            for ii = 1:obj.nel
                r = vecnorm(ec - ec(ii, :), 2, 2);
                ix = find(r < r0);
                w = 1-r(ix)/r0;
                I = [I; [ii*ones(length(ix), 1), ix, w/sw]];
            end
            Z = sparse(I(:, 1), I(:, 2), I(:, 3), obj.nel, obj.nel);
        end

        %% Misc. Function
        function plotFigs(obj)
            vM = zeros(obj.nel, 1);
            for ii = 1:obj.nel
                ix = (ii-1)*4+1:ii*4;
                vM(ii) = sqrt(obj.sigy0^2*trace(obj.sig(ix, :)*obj.P*obj.sig(ix, :)')/4);
            end
            figure;
            title("von Mises stress")
            patch(obj.ex', obj.ey', vM);
            colormap jet;
            colorbar;

            figure;
            title("Plastic Zone")
            pl = obj.ep>0;
            pl = reshape(pl, 4, obj.nel)';
            pl = any(pl, 2);
            patch(obj.ex', obj.ey', int8(pl));
            
            easyjet = [linspace(0.1, 1, 256)', ...
                       linspace(0.1, 0, 256)', ...
                       linspace(0.8, 0.1, 256)'];
            colormap(easyjet)
        end
        
        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end