classdef Solver
    properties
        edof; ex; ey
        ed; a; K
        A; t; ngp; tgp; Amax
        ndof; nel; endof; pdof; fdof 
        bcS; disp
        De; Ds; Dt; X; iX; Gam
        dDsdep; dR2dep; epst
        H; sigy0; Kinf; xi
        sigy; rtol; DP
        P; C; DeP
        R1; R1tol; N
        eps; sig; ep
        epsi; sigi; epi; Dsi; sigyi
        Z; del; p; q; ncon
        xtol; iterMax
        gam; phi; g0; g1
        sig1N
        sigTrack
        pltoel
    end

    methods
        function obj = Solver(p)  % Constructor
            % Mesh
            [~, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le);
            obj.ndof = 2*((p.lx/p.le + 1)*(p.ly/p.le + 1));
            obj.nel = round(p.lx*p.ly/p.le^2);
            obj.endof = 8;

            obj.Z = filterMatrix(obj, p.le, p.re);
            obj.p = p.p;
            obj.q = p.q;
            obj.del = p.del;
            obj.ncon = p.ncon;
            obj.xtol = p.xtol;
            obj.iterMax = p.iterMax;
            obj.g0 = zeros(obj.iterMax, 1);
            obj.g1 = zeros(obj.iterMax, 1);
            
            obj.A = p.le^2*ones(obj.nel, 1);
            obj.t = p.t;
            obj.ngp = p.ngp;
            obj.tgp = obj.nel*obj.ngp;
            obj.Amax = sum(p.Vf*obj.A);

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
            obj.xi = p.xi;

            G12 = p.E1/(2*(1+p.v12));
            obj.C = [1/p.E1, -p.v21/p.E2, -p.v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -p.v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            [obj.X, obj.iX, obj.Gam] = diagDs(obj);
            obj.Ds = repmat(obj.De, obj.tgp, 1);
            obj.Dsi = obj.Ds;
            obj.Dt = obj.Ds;

            obj.dDsdep = zeros(size(obj.Ds)); 
            obj.dR2dep = zeros(obj.tgp, 1); 
            obj.epst = zeros(obj.tgp, 1);

            obj.DeP = obj.De*obj.P;

            obj.rtol = p.rtol;
            obj.DP = p.DP;
            obj.R1tol = p.R1tol;
            obj.N = p.N;   
            obj.R1 = sparse(obj.ndof, 1);
            obj.disp = p.disp;
            obj.disp(:, 2) = obj.disp(:, 2)/obj.N;
            obj.a = zeros(obj.ndof, 1);
            obj.bcS = obj.addBC(bc, p.ly, p.le, obj.ndof);
            obj.pdof = [obj.bcS(:, 1); obj.disp(:, 1)];
            obj.fdof = setdiff((1:obj.ndof)', obj.pdof);

            obj.eps = zeros(obj.tgp, 4);
            obj.sig = zeros(obj.tgp, 4);
            obj.ep = zeros(obj.tgp, 1);
            obj.sigy = obj.sigy0*ones(obj.tgp, 1);

            obj.epsi = zeros(obj.tgp, 4);
            obj.sigi = zeros(obj.tgp, 4);
            obj.epi = zeros(obj.tgp, 1);
            obj.sigyi = obj.sigy0*ones(obj.tgp, 1);

            obj.sig1N = zeros(obj.tgp,8);
            obj.sigTrack = zeros(4,obj.N);
            obj.pltoel = zeros(obj.tgp, 2);
        end

        %% FEM
        function obj = newt(obj) % Newton-Raphson method
            for n = 1:obj.N
                % fprintf("Load step: %i \n", n);
                bc = [obj.bcS; obj.disp];
                Nr = 0;
                while norm(obj.R1(obj.fdof)) > obj.R1tol || Nr == 0
                    Nr = Nr + 1;
                    if Nr == 9
                        warning("NR converging slowly")
                    elseif Nr == 20
                        error("NR not converging")
                    end
                    obj = FEM(obj, bc, n);
                    bc(:, 2) = bc(:, 2)*0;
                    % fprintf("  Nr: %i, R1: %4.2g \n", [Nr, norm(obj.R1(obj.fdof))]);
                end
                obj.eps = obj.epsi; obj.sig = obj.sigi; obj.ep = obj.epi; obj.Ds = obj.Dsi; obj.sigy = obj.sigyi;
                obj.pltoel(:, 1) = zeros(obj.tgp, 1);
                el = 72;
                obj.sigTrack(:, n) = sqrt(obj.sigy0^2*diag(obj.sig(obj.ngp*(el-1) + 1:obj.ngp*el, :)*obj.P*obj.sig(obj.ngp*(el-1) + 1:obj.ngp*el, :)'));
                if n == 1
                    obj.sig1N(:,1:obj.ngp) = obj.sig;
                end
            end
            obj.sig1N(:,obj.ngp+1:end) = obj.sig;
        end

        function obj = FEM(obj, bc, n) % Main FEM function
            obj.K = assemK(obj, obj.Dt);

            da = obj.solvelin(obj.K, -obj.R1, bc);
            obj.a = obj.a + da;
            obj.ed = obj.a(obj.edof);

            tripf = zeros(obj.nel*obj.endof, 2);
            for el = 1:obj.nel
                fein = zeros(obj.endof, 1);
                for gp = 1:obj.ngp
                    ix = obj.ngp*(el-1) + gp;
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    [B, J] = NablaB(obj, gp, el);
                    obj.epsi(ix, [1 2 4]) = B*obj.ed(el, :)';
                    deps = (obj.epsi(ix, :) - obj.eps(ix, :))';
                    sigtr = obj.gam(el)*obj.De*deps + obj.sig(ix, :)';

                    if sqrt(obj.sigy0^2*sigtr'*obj.P*sigtr) > obj.phi(el)*obj.sigy(ix)
                            obj.pltoel(ix, 1) = n; %Ny
                        if obj.DP
                            [obj.sigi(ix, :), obj.Dt(ixM, :), obj.Dsi(ixM, :), obj.epi(ix), obj.dDsdep(ixM, :), obj.dR2dep(ix), obj.epst(ix)]...
                             = DPMat(obj, obj.epsi(ix, :)', obj.Ds(ixM, :), obj.ep(ix), obj.gam(el), obj.phi(el));
                        else
                            error("Not implemented SIMP for EPMat")
                            [obj.sigi(ix, :), obj.Dt(ixM, :), obj.sigyi(ix), Dep] = EPMat(obj, sigtr, obj.sigy(ix));
                            obj.epi(ix) = obj.epi(ix) + Dep;
                        end
                    else
                        obj.sigi(ix, :) = sigtr;
                        obj.Dt(ixM, :) = obj.gam(el)*obj.De;  
                        obj.Dsi(ixM, :) = obj.gam(el)*obj.De;
                        obj.epi(ix) = 0; obj.dDsdep(ixM, :) = 0; obj.dR2dep(ix) = 0;
                        obj.epst(ix) = (obj.gam(el)/obj.phi(el))^2*obj.epsi(ix, :)*obj.De*obj.P*obj.De*obj.epsi(ix, :)';
                            obj.pltoel(ix, 2) = obj.pltoel(ix, 1); %Ny
                    end
                    fein = fein + B'*obj.sigi(ix, [1 2 4])'*J*obj.t;
                end
                tripf((el-1)*obj.endof+1:el*obj.endof, :) = [obj.edof(el, :)', fein];
            end
            obj.R1 = sparse(tripf(:, 1), 1, tripf(:, 2), obj.ndof, 1);
        end

        function K = assemK(obj, D)
            tripK = zeros(obj.nel*obj.endof^2, 3);
            for el = 1:obj.nel
                ke = zeros(obj.endof);
                for gp = 1:obj.ngp
                    [B, J] = NablaB(obj, gp, el);
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    ke = ke + B'*D(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                end
                [rows, cols] = ndgrid(obj.edof(el, :));
                tripK((el-1)*obj.endof^2+1:el*obj.endof^2, :) = [rows(:), cols(:), ke(:)];
            end
            K = sparse(tripK(:, 1), tripK(:, 2), tripK(:, 3), obj.ndof, obj.ndof);
        end

        function [Bgp, detJ] = NablaB(obj, gp, el) % Spatial gradient of shape functions
            node = sqrt(1/3);
            gps = [-node -node; -node node; node node; node -node];
            etas = gps(gp, 1); xis = gps(gp, 2);
            Npz = 1/4*[etas-1, 1-etas, 1+etas, -1-etas;
                       xis-1, -1-xis, 1+xis, 1-xis];
            Jt = [Npz*obj.ex(el, :)', Npz*obj.ey(el, :)'];
            detJ = det(Jt);
            Npx = Jt\Npz;
            Bgp = zeros(3, 8);
            Bgp(1, 1:2:8) = Npx(1, :);
            Bgp(2, 2:2:8) = Npx(2, :);
            Npxf = flip(Npx);
            Bgp(3, :) = Npxf(:);
        end

        function a = solvelin(obj,K,f,bc) % Solve FE-equations
            nc = size(f, 2);
            if nargin == 2
                a = K\f;
            else
                s = K(obj.fdof,obj.fdof)\(f(obj.fdof,:)-K(obj.fdof,bc(:,1))*bc(:,2));
                a = zeros([obj.ndof, nc]);
                a(obj.fdof,:) = s;
                a(bc(:,1),:) = repmat(bc(:,2),1,nc);
            end
        end

        %% Material Functions
        function [sig, Dt, sigy, Dep] = EPMat(obj, sigtr, sigy)
            Dep = 0;
            sigt = sigtr'*obj.P*sigtr;
            r = sigy - obj.sigy0*sqrt(sigt);
            iter = 0;
            while norm(r) > obj.rtol || iter == 0
                iter = iter + 1;
                V = diag(1./diag(eye(4) + obj.sigy0^2*(Dep/(sigy + obj.H*Dep))*obj.Gam));
                U = obj.X*V*obj.iX;
                %dUdDep = -obj.X*V*obj.Gam*V*obj.iX*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2);
                dUdDep = -U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2);
                dsigtdU = 2*obj.P*U*(sigtr*sigtr');
                sigt = sigtr'*U*obj.P*U'*sigtr;
                drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dUdDep*dsigtdU);
                DDep = -r/drdDep;
                Dep = Dep + DDep;
                V = diag(1./diag(eye(4) + obj.sigy0^2*(Dep/(sigy + obj.H*Dep))*obj.Gam));
                U = obj.X*V*obj.iX;
                sigt = sigtr'*U*obj.P*U'*sigtr;
                r = sigy + obj.H*Dep - obj.sigy0*sqrt(sigt);           
                % fprintf("    iter: %i, r: %4.2g \n", [iter, norm(r)])
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

        function [sig, Dt, Ds, ep, dDsdep, drdep, epst] = DPMat(obj, eps, Ds, ep, gam, phi)
            epst = 1/phi^2*eps'*Ds*obj.P*Ds*eps;
            sige = (obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.xi*ep)));
            r = phi*(sige - obj.sigy0*sqrt(epst));
            iter = 0;
            while norm(r) > obj.rtol || iter == 0
                iter = iter + 1;
                if iter == 9
                    warning("Material converging slowly")
                elseif iter == 20
                    error("Material not converging")
                end
                Ds = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                detdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
                dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
                epst = 1/phi^2*(eps'*Ds*obj.P*Ds*eps);
                drdep = phi*(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
                Dep = -r/drdep;
                ep = ep + Dep;
                sige = (obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.xi*ep)));
                Ds = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                epst = 1/phi^2*(eps'*Ds*obj.P*Ds*eps);
                r = phi*(sige - obj.sigy0*sqrt(epst));
                % fprintf("    iter: %i, r: %4.2g \n", [iter, norm(r)])
            end
            sig = Ds*eps;
            detdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
            dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
            drdep =  phi*(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            drdeps = 1/phi*(-obj.sigy0/sqrt(epst)*Ds*obj.P*Ds*eps);
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
            I = zeros(obj.nel*(2*re)^2, 3);
            [x, y] = meshgrid(-re:re, -re:re);
            weights = max(0, 1 - sqrt(x.^2 + y.^2)/re);
            sw = sum(weights(:));
            r0 = le*re;
            i = 0;
            for ii = 1:obj.nel
                r = vecnorm(ec - ec(ii, :), 2, 2);
                ix = find(r < r0);
                ixn = length(ix);
                w = 1-r(ix)/r0;
                I(i+1:i+ixn, :) = [ii*ones(ixn, 1), ix, w/sw];
                i = i + ixn;
            end
            Z = sparse(I(1:i, 1), I(1:i, 2), I(1:i, 3), obj.nel, obj.nel);
        end

        function obj = optimizer(obj, x)
            a0 = 1; a1 = zeros(obj.ncon,1); c = 1000*ones(obj.ncon,1); d = ones(obj.ncon,1);
            xold1 = []; xold2 = []; low = []; upp = [];
            dx = 1;
            iter = 0;
            while dx > obj.xtol
                iter = iter + 1;
                if iter == obj.iterMax + 1
                    fprintf("\n\nMax iteration count reached")
                    break
                end

                obj = initOpt(obj, x);
                obj = newt(obj);
                [obj.g0(iter), dg0, obj.g1(iter), dg1] = funcEval(obj, x);
                if iter == 1
                    s = abs(100/obj.g0(1));
                end
                [xnew,~,~,~,~,~,~,~,~,low,upp] = mmasub(obj.ncon, obj.nel, iter, x, zeros(obj.nel, 1), ones(obj.nel, 1), ...
                                                            xold1, xold2, s*obj.g0(iter), s*dg0, obj.g1(iter), dg1, low, upp, a0, a1, c, d);
                xold2 = xold1;
                xold1 = x;
  
                dx = norm((xnew - x)/obj.nel);
                % x = obj.Z*xnew;
                x = xnew;

                plotFigs(obj, x, 0, 1);
                fprintf("Opt iter: %i\n", iter)
                fprintf("  g0: %.2g, g1: %.2g, dx: %.2g\n", [obj.g0(iter), obj.g1(iter), dx])
            end
            obj.g0 = obj.g0(1:iter);
            obj.g1 = obj.g1(1:iter);
            plotFigs(obj, x, iter, 0);
        end

        function [g0, dg0, g1, dg1] = funcEval(obj, x)
            dgt0dx = zeros(1, obj.nel);
            dR1dx = zeros(obj.ndof, obj.nel);
            dR2dxe = zeros(obj.ngp, 1);
            dR2dx = zeros(obj.tgp, obj.nel);
            dgt0dep = zeros(obj.tgp, 1);
            dR1dep = zeros(obj.ndof, obj.tgp);
            ap = zeros(obj.ndof, 1);
            ap(obj.pdof) = obj.a(obj.pdof);

            for el = 1:obj.nel
                dgam = obj.p*(1-obj.del)*x(el)^(obj.p-1);
                dphi = obj.q*(1-obj.del)*x(el)^(obj.q-1);
                th = (dgam*obj.phi(el)-dphi*obj.gam(el))/obj.phi(el)^2;
                Kte = zeros(obj.endof);
                eix = obj.edof(el, :);
                for gp = 1:obj.ngp
                    [B, J] = NablaB(obj, gp, el);
                    ix = obj.ngp*(el-1) + gp;
                    ixM =  4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    k0 = obj.ep(ix)*obj.sigy0^2/(obj.sigy0+obj.H*obj.ep(ix));
                    % V = inv(obj.De + obj.gam(el)/obj.phi(el)*k0*obj.De*obj.P*obj.De);
                    % dDsdx = dgam*obj.De*V*obj.De - obj.gam(el)*obj.De*V*(th*k0*obj.De*obj.P*obj.De)*V*obj.De;
                    dDsdx = (dgam*obj.Ds(ixM, :) - th*k0*obj.Ds(ixM, :)*obj.P*obj.Ds(ixM, :))/obj.gam(el);
                    Kte = Kte + B'*dDsdx(([1 2 4]),[1 2 4])*B*J*obj.t;

                    depstdx = obj.eps(ix, :)*(dDsdx*obj.P*obj.Ds(ixM, :)...
                              - obj.Ds(ixM, :)*obj.P*2*dphi/obj.phi(el)*obj.Ds(ixM, :)...
                              + obj.Ds(ixM, :)*obj.P*dDsdx)*obj.eps(ix, :)'/obj.phi(el)^2;
                    dR2dxe(gp) = dphi*(obj.sigy0 + obj.H*obj.ep(ix))...
                                 - dphi*obj.sigy0*sqrt(obj.epst(ix)) - obj.phi(el)*obj.sigy0/2/sqrt(obj.epst(ix))*depstdx;

                    Kh = B'*obj.dDsdep(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                    dR1depe = Kh*obj.a(eix);
                    dgt0dep(ix) = -ap(eix)'*dR1depe;
                    dR1dep(eix, ix) = dR1depe;
                end
                dR1dxe = Kte*obj.a(eix);
                dgt0dx(el) = -ap(eix)'*dR1dxe;
                dR1dx(eix, el) = dR1dxe;
                dR2dx(ix-3:ix, el) = dR2dxe;
            end

            dgt0du = -obj.a(obj.pdof)'*obj.K(obj.pdof, obj.fdof);

            lamt = -dgt0du/obj.K(obj.fdof, obj.fdof);
            idR2dep = diag(1./obj.dR2dep);
            idR2dep(isinf(idR2dep)) = 0;
            mut = -dgt0dep'*idR2dep - lamt*dR1dep(obj.fdof, :)*idR2dep;

            g0 = -obj.a(obj.pdof)'*obj.R1(obj.pdof);
            dg0 = obj.Z'*(dgt0dx + lamt*dR1dx(obj.fdof, :) + mut*dR2dx)';
            % dg0 = (dgt0dx + lamt*dR1dx(obj.fdof, :) + mut*dR2dx)';

            g1 = x'*obj.A/obj.Amax - 1;
            % dg1 = (obj.Z'*obj.A/obj.Amax)';
            dg1 = (obj.A/obj.Amax)';
        end

        function obj = initOpt(obj, x)
            obj.eps = zeros(obj.tgp, 4);
            obj.sig = zeros(obj.tgp, 4);
            obj.ep = zeros(obj.tgp, 1);
            obj.epsi = zeros(obj.tgp, 4);
            obj.sigi = zeros(obj.tgp, 4);
            obj.epi = zeros(obj.tgp, 1);

            obj.a = zeros(obj.ndof, 1);
            obj.R1 = sparse(obj.ndof, 1);

            
            obj.gam = obj.del + (1-obj.del)*x.^obj.p;
            obj.phi = obj.del + (1-obj.del)*x.^obj.q;
            gam4 = repelem(obj.gam, 4*obj.ngp);
            obj.Ds = gam4.*repmat(obj.De, obj.tgp, 1);
            obj.Dsi = obj.Ds;
            obj.Dt = obj.Ds;

            obj.dDsdep = zeros(size(obj.Ds)); 
            obj.dR2dep = zeros(obj.tgp, 1);
        end

        

        %% Misc. Function
        function plotFigs(obj, x, iter, flag)
            if flag
                clf;
                colormap(flipud(gray(256)));
                patch(obj.ex', obj.ey', x);
                colorbar;
                axis equal;
                drawnow;
            else
                vM = zeros(obj.nel, 1);
                for ii = 1:obj.nel
                    ix = (ii-1)*4+1:ii*4;
                    vM(ii) = sqrt(obj.sigy0^2*trace(obj.sig(ix, :)*obj.P*obj.sig(ix, :)')/4);
                end
                
                cosT = zeros(obj.nel,1);
                for ii = 1:obj.nel
                    ix = (ii-1)*4+1:ii*4;
                    cosT(ii) = trace(obj.sig1N(ix,1:obj.ngp)*obj.sig1N(ix,obj.ngp+1:end)')/(vecnorm(obj.sig1N(ix,1:obj.ngp),2,2)'*vecnorm(obj.sig1N(ix,obj.ngp+1:end),2,2));
                end

                
                figure;
                title("cos(\theta) field for compliance maximization");
                patch(obj.ex', obj.ey', cosT);
                colormap jet;
                colorbar;

                figure;
                title("von Mises stress");
                patch(obj.ex', obj.ey', vM);
                colormap jet;
                colorbar;
    
                figure;
                title("Plastic Zone");
                pl = obj.ep>0;
                pl = reshape(pl, 4, obj.nel)';
                pl = any(pl, 2);
                patch(obj.ex', obj.ey', int8(pl));
                
                softjet = [linspace(0.1, 1, 256)', ...
                           linspace(0.1, 0, 256)', ...
                           linspace(0.8, 0.1, 256)'];
                colormap(softjet);

                figure;
                plot(1:iter, obj.g0, LineWidth=2)
                title("Convergance Plot g0");

                figure;
                plot(1:iter, obj.g1, LineWidth=2)
                title("Convergance Plot g1")

            end
        end
        
        function bc = addBC(~, bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end