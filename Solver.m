classdef Solver
    properties
        edof; ex; ey; axi
        ed; a; K
        A; t; ngp; tgp; Amax
        ndof; nel; endof; pdof; fdof 
        bcS; disp; fixDens
        Bgp; detJ; rows; cols
        De; Ds; Dt; X; iX; Gam
        dDsdep; dR2dep; epst
        H; sigy0; Kinf; xi
        sigy; rtol; PT
        P; C; DeP
        R1; R1tol; N
        eps; sig; ep
        epsi; sigi; epi; Dsi; sigyi
        Z; filtOn; eta; beta
        del; dels; p; q; ncon
        rampB; rampPQ
        xtol; iterMax
        gam; phi; g0; g1
        sig1N
        saveName; prints
    end

    methods
        function obj = Solver(p)  % Constructor
            [~, ~, ~, obj.edof, obj.ex, obj.ey, bc] = designDomain(p.lx, p.ly, p.le, p.wx, p.wy);
            obj.nel = length(obj.edof(:,1));
            obj.ndof = max(max(obj.edof(:,2:end)));
            obj.endof = 8;
            obj.N = p.N;
            obj.axi = [0 p.lx 0 p.ly];

            if p.loadcase == 1 % Vertical load on beam
                obj.disp(:, 1) = 2:2:10;
                obj.disp(:, 2) = p.disp*ones(5,1);
                obj.bcS = obj.addBC(bc, p.ly, p.le, obj.ndof);
            elseif p.loadcase == 2 % Bending beam
                obj.disp(:, 1) = bc(:,1);
                obj.disp(:, 2) = linspace(-p.disp, p.disp, length(bc));
                obj.bcS = obj.addBC([], p.ly, p.le, obj.ndof);
            elseif p.loadcase == 3 % L-beam
                obj.disp(:,1) =  bc(bc(:,2) == 1, 1);
                obj.disp = [obj.disp, p.disp*ones(size(obj.disp))];
                obj.bcS = bc(bc(:,2) == 0, 1);
                obj.bcS = [obj.bcS, zeros(size(obj.bcS))];
            else
                error("Load case doesn't exist");
            end
            obj.disp(:, 2) = obj.disp(:, 2)/obj.N;
            % obj.fixDens = find(any(ismember(obj.edof,obj.disp(:,1)),2));
            
            obj.saveName = p.saveName;
            obj.prints = p.print;
            obj.filtOn = p.filtOn;
            obj.Z = filterMatrix(obj, p.le, p.re);
            obj.p = p.p;
            obj.q = p.q;
            obj.eta = p.eta;    
            obj.beta = p.beta;
            obj.rampB = p.rampB;
            obj.rampPQ = p.rampPQ;
            
            obj.del = p.del;
            obj.dels = p.dels;
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
            obj.Bgp = zeros(3*obj.ngp, obj.endof);
            obj.detJ = zeros(obj.ngp, 1);

            for gp = 1:obj.ngp
                [obj.Bgp(3*(gp-1)+1:3*gp, :), obj.detJ(gp)] = nablaB(obj, gp, 1);
            end
            obj.rows = zeros(obj.endof^2, obj.nel);
            obj.cols = zeros(obj.endof^2, obj.nel);
            for el = 1:obj.nel
                [rows, cols] = ndgrid(obj.edof(el, :));
                obj.rows(:, el) = rows(:);
                obj.cols(:, el) = cols(:);
            end

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
            v21 = p.v12 * (p.E2 / p.E1);  % From reciprocity
            v31 = p.v13 * (p.E3 / p.E1);
            v32 = p.v23 * (p.E3 / p.E2);
            obj.C = [1/p.E1, -v21/p.E2, -v31/p.E3, 0;
                    -p.v12/p.E1, 1/p.E2, -v32/p.E3, 0;
                    -p.v13/p.E1, -p.v23/p.E2, 1/p.E3, 0;
                     0, 0, 0, 1/G12];
            obj.De = inv(obj.C);
            [obj.X, obj.iX, obj.Gam] = diagDe(obj);
            obj.Ds = repmat(obj.De, obj.tgp, 1);
            obj.Dsi = obj.Ds;
            obj.Dt = obj.Ds;

            obj.dDsdep = zeros(size(obj.Ds)); 
            obj.dR2dep = zeros(obj.tgp, 1); 
            obj.epst = zeros(obj.tgp, 1);

            obj.DeP = obj.De*obj.P;

            obj.rtol = p.rtol;
            obj.PT = p.PT;
            obj.R1tol = p.R1tol;
            obj.R1 = sparse(obj.ndof, 1);
         
            obj.a = zeros(obj.ndof, 1);

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
        end

        %% Optimization
        function [obj, x] = opt(obj, x) % Optimizer with MMA
            a0 = 1; a1 = zeros(obj.ncon,1); c = 1000*ones(obj.ncon,1); d = ones(obj.ncon,1);
            xold1 = []; xold2 = []; low = []; upp = [];
            x(obj.fixDens) = 1;
            dx = 1;
            iter = 0;
            while dx > obj.xtol
                iter = iter + 1;
                if iter == obj.iterMax + 1
                    iter = obj.iterMax;
                    fprintf("\n\nMax iteration count reached\n")
                    break
                elseif mod(iter, 10) == 0
                    if obj.rampB == 1 && obj.beta < 10
                        obj.beta = obj.beta*1.1;
                    elseif obj.rampB == 2 && obj.beta < 10
                        obj.beta = obj.beta + 1;
                    end
                    if obj.rampPQ && obj.p < 3
                        obj.p = obj.p + 0.1;
                        obj.q = obj.q + 0.1;
                    end
                end
                obj = init(obj, x);
                obj = newt(obj);
                [obj.g0(iter), dg0, obj.g1(iter), dg1] = funcEval(obj, x);
                if iter == 1
                    s = abs(100/obj.g0(1));
                end
                [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(obj.ncon, obj.nel, iter, x, zeros(obj.nel, 1), ones(obj.nel, 1), ...
                                                            xold1, xold2, s*obj.g0(iter), s*dg0, obj.g1(iter), dg1, low, upp, a0, a1, c, d);
                xold2 = xold1;
                xold1 = x;
  
                dx = norm((xmma - x)/obj.nel);
                x = xmma;

                plotFigs(obj, x, 1);
                fprintf("Opt iter: %i\n", iter)
                fprintf("  g0: %.2g, g1: %.2g, dx: %.2g\n", [obj.g0(iter), obj.g1(iter), dx])
            end
            obj.g0 = obj.g0(1:iter);
            obj.g1 = obj.g1(1:iter);
            plotFigs(obj, x, 0);
        end

        function [g0, dg0, g1, dg1] = funcEval(obj, x)
            dgt0dx = zeros(1, obj.nel);
            dR1dx = zeros(obj.endof*obj.nel, 1);  
            dR2dx = zeros(obj.tgp, 1);  
            dgt0dep = zeros(obj.tgp, 1);
            dR1dep = zeros(obj.endof*obj.tgp, 1);
            ap = sparse(obj.pdof, 1, obj.a(obj.pdof), obj.ndof, 1);

            x = obj.Z*x;
            [x, dxH] = he(obj,x);
            for el = 1:obj.nel
                dgam = obj.p*(1-obj.del)*x(el)^(obj.p-1);
                dphi = obj.q*(1-obj.dels)*x(el)^(obj.q-1);
                th = (dgam*obj.phi(el)-dphi*obj.gam(el))/obj.phi(el)^2;
                Kte = zeros(obj.endof);
                eix = obj.edof(el, :);
                for gp = 1:obj.ngp
                    B = obj.Bgp(3*(gp-1)+1:3*gp, :); J = obj.detJ(gp); % [B, J] = NablaB(obj, gp, el);
                    ix = obj.ngp*(el-1) + gp;
                    ixM =  4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    k0 = obj.ep(ix)*obj.sigy0^2/(obj.sigy0 + obj.H*obj.ep(ix) + obj.Kinf*(1-exp(-obj.xi*obj.ep(ix))));
                    dDsdx = (dgam*obj.Ds(ixM, :) - th*k0*obj.Ds(ixM, :)*obj.P*obj.Ds(ixM, :))/obj.gam(el);
                    Kte = Kte + B'*dDsdx(([1 2 4]),[1 2 4])*B*J*obj.t;

                    depstdx = obj.eps(ix, :)*(dDsdx*obj.P*obj.Ds(ixM, :)...
                              - obj.Ds(ixM, :)*obj.P*2*dphi/obj.phi(el)*obj.Ds(ixM, :)...
                              + obj.Ds(ixM, :)*obj.P*dDsdx)*obj.eps(ix, :)'/obj.phi(el)^2;
                    dR2dx(ix) = dphi*(obj.sigy0 + obj.H*obj.ep(ix) + obj.Kinf*(1-exp(-obj.xi*obj.ep(ix))))...
                                 - dphi*obj.sigy0*sqrt(obj.epst(ix)) - obj.phi(el)*obj.sigy0/2/sqrt(obj.epst(ix))*depstdx;

                    Kh = B'*obj.dDsdep(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                    dR1depe = Kh*obj.a(eix);
                    dgt0dep(ix) = -ap(eix)'*dR1depe;
                    dR1dep(8*(ix-1)+1:8*ix) = dR1depe; 
                end
                dR1dxe = Kte*obj.a(eix);
                dgt0dx(el) = -ap(eix)'*dR1dxe;
                dR1dx(8*(el-1)+1:8*el) = dR1dxe; 
            end
            dR1dx = sparse(reshape(obj.edof', [], 1), repelem((1:obj.nel)', obj.endof), dR1dx, obj.ndof, obj.nel);
            dR2dx = sparse((1:obj.tgp)', repelem((1:obj.nel)', obj.ngp), dR2dx, obj.tgp, obj.nel);
            dR1dep = sparse(reshape(repelem(obj.edof', 1, obj.ngp), [], 1), repelem((1:obj.tgp)', obj.endof), dR1dep, obj.ndof, obj.tgp);

            dgt0da = -obj.a(obj.pdof)'*obj.K(obj.pdof, obj.fdof);

            pgp = find(obj.ep);
            lamt = -dgt0da/obj.K(obj.fdof, obj.fdof);
            idR2dep = diag(1./obj.dR2dep(pgp));
            mut = -dgt0dep(pgp)'*idR2dep - lamt*dR1dep(obj.fdof, pgp)*idR2dep;

            g0 = -obj.a(obj.pdof)'*obj.R1(obj.pdof);
            dg0 = dxH'.*obj.Z'*(dgt0dx + lamt*dR1dx(obj.fdof, :) + mut*dR2dx(pgp, :))';
            dg0(obj.fixDens) = 0;

            g1 = x'*obj.A/obj.Amax - 1;
            dg1 = (dxH'.*obj.Z'*obj.A/obj.Amax)';
            dg1(obj.fixDens) = 0;
        end

        function Z = filterMatrix(obj, le, re)
            if obj.filtOn
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
            else
                Z = speye(obj.nel);
            end
        end

        function [x, dxH] = he(obj, x)
            if obj.filtOn
                dxH = obj.beta*(1-tanh(obj.beta*(x-obj.eta)).^2)/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));
                x = (tanh(obj.beta*obj.eta) + tanh(obj.beta*(x-obj.eta)))/(tanh(obj.beta*obj.eta) + tanh(obj.beta*(1-obj.eta)));
            else
                dxH = ones(obj.nel, 1);
            end
        end

        %% FEM
        function obj = newt(obj) % Newton-Raphson method
            tdisp = abs(obj.disp(1, 2)*obj.N); cdisp = obj.disp; idisp = 0; pn = obj.N;
            rst = 0; f = 1;
            n = 1;
            while abs(tdisp-idisp) > 1e-8
                if obj.prints(1)
                    fprintf("Load step: %i \n", n);
                end
                bc = [obj.bcS; cdisp];
                Nr = 0;
                while norm(obj.R1(obj.fdof)) > obj.R1tol || Nr == 0
                    Nr = Nr + 1;
                    if Nr == 9
                        rst = rst + 1;
                        if rst == 7
                            f = 30;
                        elseif rst == 20
                            error("Too many restarts");
                        end
                        pn = nextprime(pn+f);
                        cdisp(:, 2) = obj.disp(:, 2)*obj.N/pn;
                        obj = init(obj);
                        idisp = 0; n = 1;
                        warning("Restarting Nr %i", rst);
                        break;
                    end
                    obj = FEM(obj, bc);
                    bc(:, 2) = bc(:, 2)*0;
                    if obj.prints(2)
                        fprintf("  Nr: %i, R1: %4.2g \n", [Nr, norm(obj.R1(obj.fdof))]);
                    end
                end
                if Nr == 9
                    continue;
                end
                idisp = idisp + abs(cdisp(1, 2));
                obj.eps = obj.epsi; obj.sig = obj.sigi; obj.ep = obj.epi; obj.Ds = obj.Dsi; obj.sigy = obj.sigyi;
                if n == 1
                    obj.sig1N(:,1:obj.ngp) = obj.sig;
                end
                n = n + 1;
            end
            obj.sig1N(:,obj.ngp+1:end) = obj.sig;
        end

        function obj = FEM(obj, bc) % Main FEM function
            obj.K = assemK(obj, obj.Dt);
            da = obj.solveLin(obj.K, -obj.R1, bc);
            obj.a = obj.a + da;
            obj.ed = obj.a(obj.edof);

            tripf = zeros(obj.nel*obj.endof, 2);
            for el = 1:obj.nel
                fein = zeros(obj.endof, 1);
                for gp = 1:obj.ngp
                    ix = obj.ngp*(el-1) + gp;
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    B = obj.Bgp(3*(gp-1)+1:3*gp, :); J = obj.detJ(gp); % [B, J] = NablaB(obj, gp, el);
                    obj.epsi(ix, [1 2 4]) = B*obj.ed(el, :)';
                    deps = (obj.epsi(ix, :) - obj.eps(ix, :))';
                    sigtr = obj.gam(el)*obj.De*deps + obj.sig(ix, :)';

                    if sqrt(obj.sigy0^2*sigtr'*obj.P*sigtr) > obj.phi(el)*obj.sigy(ix)
                        if obj.PT
                            [obj.sigi(ix, :), obj.Dt(ixM, :), obj.Dsi(ixM, :), obj.epi(ix), obj.dDsdep(ixM, :), obj.dR2dep(ix), obj.epst(ix)]...
                             = DPMat(obj, obj.epsi(ix, :)', obj.Ds(ixM, :), obj.ep(ix), obj.gam(el), obj.phi(el));
                        else
                            [obj.sigi(ix, :), obj.Dt(ixM, :), obj.sigyi(ix), Dep] = IPMat(obj, sigtr, obj.sigy(ix), obj.gam(el), obj.phi(el));
                            obj.epi(ix) = obj.epi(ix) + Dep;
                        end
                    else
                        obj.sigi(ix, :) = sigtr;
                        obj.Dt(ixM, :) = obj.gam(el)*obj.De;  
                        obj.Dsi(ixM, :) = obj.gam(el)*obj.De;
                        obj.epi(ix) = 0; obj.dDsdep(ixM, :) = 0; obj.dR2dep(ix) = 0;
                        obj.epst(ix) = (obj.gam(el)/obj.phi(el))^2*obj.epsi(ix, :)*obj.De*obj.P*obj.De*obj.epsi(ix, :)';
                    end
                    fein = fein + B'*obj.sigi(ix, [1 2 4])'*J*obj.t;
                end
                tripf((el-1)*obj.endof+1:el*obj.endof, :) = [obj.edof(el, :)', fein];
            end
            obj.R1 = sparse(tripf(:, 1), 1, tripf(:, 2), obj.ndof, 1);
        end

        function K = assemK(obj, D) % Assemble stiffness matrix
            tripK = zeros(obj.nel*obj.endof^2, 3);
            for el = 1:obj.nel
                ke = zeros(obj.endof);
                for gp = 1:obj.ngp
                    B = obj.Bgp(3*(gp-1)+1:3*gp, :); J = obj.detJ(gp); % [B, J] = nablaB(obj, gp, el);
                    ixM = 4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    ke = ke + B'*D(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                end
                tripK((el-1)*obj.endof^2+1:el*obj.endof^2, :) = [obj.rows(:, el), obj.cols(:, el), ke(:)];
            end
            K = sparse(tripK(:, 1), tripK(:, 2), tripK(:, 3), obj.ndof, obj.ndof);
        end

        function [Bgp, detJ] = nablaB(obj, gp, el) % Spatial gradient of shape functions
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

        function a = solveLin(obj,K,f,bc) % Solve FE-equations
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
        function [sig, Dt, Ds, ep, dDsdep, drdep, epst] = DPMat(obj, eps, Ds, ep, gam, phi) % Deformation Plasticity
            epst = 1/phi^2*eps'*Ds*obj.P*Ds*eps;
            sige = (obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.xi*ep)));
            r = (sige - obj.sigy0*sqrt(epst));
            iter = 0;
            while norm(r) > obj.rtol || iter == 0
                iter = iter + 1;
                if iter == 9
                    warning("Material converging slowly")
                elseif iter == 15
                    error("Material not converging")
                end
                detdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
                dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
                drdep = (obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
                Dep = -r/drdep;
                ep = ep + Dep;
                sige = (obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.xi*ep)));
                Ds = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
                epst = 1/phi^2*(eps'*Ds*obj.P*Ds*eps);
                r = (sige - obj.sigy0*sqrt(epst));
                if obj.prints(3)
                    fprintf("    iter: %i, r: %4.2g \n", [iter, norm(r)])
                end
            end
            sig = Ds*eps;
            detdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
            dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
            drdep =  phi*(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            drdeps = 1/phi*(-obj.sigy0/sqrt(epst)*Ds*obj.P*Ds*eps);
            depdeps = -drdeps/drdep;
            Dt = Ds + dDsdep*eps*depdeps';
        end

        function [sig, Dt, sigy, Dep] = IPMat(obj, sigtr, sigy, gam, phi) % Incremental Plasticity
            Dep = 0;
            U = eye(4);
            sigt = 1/phi^2*(sigtr'*obj.P*sigtr);
            r = sigy - obj.sigy0*sqrt(sigt);
            iter = 0;
            while norm(r) > obj.rtol || iter == 0
                iter = iter + 1;
                if iter == 9
                    warning("Material converging slowly")
                elseif iter == 15
                    error("Material not converging")
                end
                dsigtdU = 1/phi^2*(2*obj.P*U*(sigtr*sigtr'));
                dUdDep = gam/phi*(-U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2));
                drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dsigtdU'*dUdDep);
                DDep = -r/drdDep;
                Dep = Dep + DDep;
                U = obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/(sigy + obj.H*Dep)*Dep*obj.Gam))*obj.iX;
                sigt = 1/phi^2*(sigtr'*U'*obj.P*U*sigtr);
                r = sigy + obj.H*Dep - obj.sigy0*sqrt(sigt);
                if obj.prints(3)
                    fprintf("    iter: %i, r: %4.2g \n", [iter, norm(r)])
                end
            end
            sig = U*sigtr;
            sigy = sigy + obj.H*Dep;
            dUdDep = gam/phi*(-U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2));
            dsigtdU = 1/phi^2*(2*obj.P*U*(sigtr*sigtr'));
            drdDep = phi*(obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dsigtdU'*dUdDep));
            drdDeps = gam/phi*(-obj.sigy0/sqrt(sigt)*obj.De*U'*obj.P*U*sigtr);
            dDepdDeps = -drdDeps/drdDep;
            Dt = gam*U*obj.De + dUdDep*sigtr*dDepdDeps';
        end

        function [X, iX, Gam] = diagDe(obj)
            [Q, L] = svd(obj.De);
            sL = sqrt(L);
            B = sL*Q'*obj.P*Q*sL;
            [R, Gam] = svd(B);
            X = Q*sL*R;
            iX = inv(X);
        end

        %% Misc. Function
        function obj = init(obj, x)
            if nargin == 2
                x = he(obj, obj.Z*x);
                obj.gam = obj.del + (1-obj.del)*x.^obj.p;
                obj.phi = obj.dels + (1-obj.dels)*x.^obj.q;
            end
            gam4 = repelem(obj.gam, 4*obj.ngp);
            obj.Ds = gam4.*repmat(obj.De, obj.tgp, 1);
            obj.Dsi = obj.Ds;
            obj.Dt = obj.Ds;

            obj.eps = zeros(obj.tgp, 4);
            obj.sig = zeros(obj.tgp, 4);
            obj.ep = zeros(obj.tgp, 1);
            obj.epsi = zeros(obj.tgp, 4);
            obj.sigi = zeros(obj.tgp, 4);
            obj.epi = zeros(obj.tgp, 1);

            obj.a = zeros(obj.ndof, 1);
            obj.R1 = sparse(obj.ndof, 1);

            obj.dDsdep = zeros(size(obj.Ds)); 
            obj.dR2dep = zeros(obj.tgp, 1);
        end

        function plotFigs(obj, x, flag)
            if flag
                clf;
                colormap(flipud(gray(256)));
                patch(obj.ex', obj.ey', x, ...
                    'EdgeColor', 'none');
                axis equal;
                axis(obj.axi);
                axis off;
                drawnow;
            else
                cosT = zeros(obj.nel,1);
                Hs = zeros(obj.nel, 1);
                for el = 1:obj.nel
                    ix = obj.ngp*(el-1)+1:el*obj.ngp;
                    cosT(el) = trace(obj.sig1N(ix,1:obj.ngp)*obj.sig1N(ix,obj.ngp+1:end)')/(vecnorm(obj.sig1N(ix,1:obj.ngp),2,2)'*vecnorm(obj.sig1N(ix,obj.ngp+1:end),2,2));
                    Hs(el) = sqrt(obj.sigy0^2*trace(obj.sig(ix, :)*obj.P*obj.sig(ix, :)')/obj.ngp);
                end
                cosT(x<0.5) = 0;
                Hc = Hs./(obj.phi*obj.sigy0);

                figure;
                tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
                nexttile;
                title("cos($\theta$) field for compliance maximization", 'Interpreter', 'latex');
                patch(obj.ex', obj.ey', cosT, ...
                      'EdgeColor', 'none');
                colormap(gca, 'jet'); 
                clim([0 1]);
                colorbar('FontSize', 16);
                axis equal off;
                axis(obj.axi);

                nexttile;
                title("Effective Hill Stress", 'Interpreter', 'latex');
                patch(obj.ex', obj.ey', Hs*1e-6, ...
                      'EdgeColor', 'none');
                colormap(gca, 'jet');
                colorbar('FontSize', 16);
                axis equal off;
                axis(obj.axi);

                figure;
                title("Plastic Zone", 'Interpreter', 'latex');
                pl = obj.ep>0;
                pl = reshape(pl, 4, obj.nel)';
                pl = any(pl, 2);
                patch(obj.ex(pl, :)', obj.ey(pl, :)', 'red', ...
                      'EdgeColor', 'none', ...
                      'DisplayName', 'Plasticity');
                patch(obj.ex(~pl, :)', obj.ey(~pl, :)', 'blue', ...
                      'EdgeColor', 'none', ...
                      'DisplayName', 'Elasticity');
                patch(obj.ex(x<0.5, :)', obj.ey(x<0.5, :)', 'white', ...
                      'EdgeColor', 'none', ...
                      'HandleVisibility', 'off');
                legend('Location', 'south', 'FontSize', 12)
                axis equal off;
                axis(obj.axi);

                figure;
                title("Hill Criterion", 'Interpreter', 'latex');
                patch(obj.ex', obj.ey', Hc, ...
                      'EdgeColor', 'none');
                colormap(gca, 'jet');
                clim([0 1.5]);
                colorbar('FontSize', 16);
                axis equal off;
                axis(obj.axi);

                iter = length(obj.g0);
                if iter < 10
                    return
                else
                    figure;
                    yyaxis left;
                    plot((10:iter)', -obj.g0(10:end), 'r:', 'LineWidth', 2);
                    ylabel('g0 objective');
                    ylim([0 round(max(-obj.g0(10:end))+100, -2)]);
                    ax = gca;
                    ax.YColor = 'k';

                    yyaxis right;
                    plot((10:iter)', obj.g1(10:end), 'k', 'LineWidth', 2);
                    ylabel('g1 constraint');
                    ylim([-0.5 0.5]);
                    ax = gca;
                    ax.YColor = 'k';

                    xlabel('Iteration');
                    legend('Stiffness', 'Volume');            
                    title("Convergence Plot");
                    grid on;
                end
            end
        end

        function saveData(obj, x, params, path)
            if obj.saveName == ""
                return
            else
                if nargin < 4
                    path = "";
                end
                dataPath = fullfile(path, sprintf("%s.mat", obj.saveName));
                iter = 0;
                while isfile(dataPath)
                    iter = iter + 1;
                    dataPath = fullfile(path, sprintf("%s_copy%i.mat", obj.saveName, iter));
                end
                val = obj.assignVar(obj, struct());             
                save(dataPath, "x", "params", "val");
            end
        end
    end

    methods (Static)
        function out = assignVar(in, out)
            fields = {'ex', 'ey', 'ngp', 'nel', 'sig', 'P', 'sigy0', 'sig1N', 'ep', 'g0', 'g1'};
            for i = 1:numel(fields)
                out.(fields{i}) = in.(fields{i});
            end
        end

        function bc = addBC(bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end
    end
end