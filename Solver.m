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
        Z; Zp; pc; filtOn; eta; beta
        stressCon; pnm; mmaVals; zeroGrad
        sigm; cp; ca
        del; dels; p; q; ncon
        rampB; rampPQ
        dx; xtol; ftol; iterMax
        gam; phi; g0; gc
        sig1N; logs
        saveName; prints; plots; design
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
                obj.bcS = bc(bc(:,2) == 0, :);
            else
                error("Load case doesn't exist");
            end
            obj.disp(:, 2) = obj.disp(:, 2)/obj.N;
            %obj.fixDens = find(any(ismember(obj.edof,obj.disp(:,1)),2));

            obj.saveName = p.saveName;
            if isfield(p, 'prints'); obj.prints = p.prints; else; obj.prints = [0,0,0]; end
            if isfield(p, 'plots'); obj.plots = p.plots; else; obj.plots = 1; end
            obj.filtOn = p.filtOn;
            obj.pc = obj.padding(p.lx, p.ly, p.le, p.wx, p.wy, p.re, p.loadcase, p.pad);
            [obj.Z, obj.Zp] = filterMatrix(obj, p.le, p.re, obj.pc);
            obj.p = p.p;
            obj.q = p.q;
            obj.eta = p.eta;
            obj.beta = p.beta;
            obj.rampB = p.rampB;
            obj.rampPQ = p.rampPQ;

            obj.stressCon = p.stressCon;
            obj.pnm = p.pnm;
            obj.cp = 1;
            obj.ca = 1;
            if isfield(p, 'mmaEnd'); obj.mmaVals = [p.mmaEnd p.mma]; else; obj.mmaVals = [0, 0.5, 10, 0.01, 0.5, 10, 0.01]; end

            obj.del = p.del;
            obj.dels = p.dels;
            obj.ncon = 1 + obj.stressCon;
            obj.xtol = p.xtol;
            if isfield(p, 'ftol'); obj.ftol = p.ftol; else; obj.ftol = 0.1; end
            obj.iterMax = p.iterMax;
            obj.dx = zeros(obj.iterMax + 1, 1);
            obj.g0 = zeros(obj.iterMax + 1, 1);
            obj.gc = zeros(obj.iterMax + 1, 1 + obj.stressCon);
            
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
            obj.sigm = p.sigc*obj.sigy0*ones(obj.nel, 1);
            obj.sigm(obj.ex(:, 1) >= p.lx - p.le*(p.stressFree + 1e-5)) = p.sigc*obj.sigy0*1e15;
            if isfield(p, 'zeroGrad'); obj.zeroGrad = p.zeroGrad & (obj.ex(:, 1) >= p.lx - p.le*(p.stressFree + 1e-5)); else; obj.zeroGrad = false; end

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
            sigynel = obj.sigy0*ones(obj.nel, 1);
            sigynel(obj.ex(:, 1) >= p.lx - p.le*(p.plasticFree + 1e-5)) = obj.sigy0*1e15;
            obj.sigy = repelem(sigynel, 4);
            
            obj.epsi = obj.eps;
            obj.sigi = obj.sig;
            obj.epi = obj.ep;
            obj.sigyi = obj.sigy;

            obj.sig1N = zeros(obj.tgp, 8);
            obj.logs = struct();
        end

        %% Optimization
        function [obj, x] = opt(obj, x) % Optimizer with MMA
            a0 = 1; a1 = zeros(obj.ncon,1); c = 1000*ones(obj.ncon,1); d = ones(obj.ncon,1);
            xold1 = []; xold2 = []; low = []; upp = [];
            x(obj.fixDens) = 1;
            iter = 0;
            criteria = true;
            simp_iter = inf;

            while criteria
                iter = iter + 1;
                if iter == obj.iterMax + 1
                    iter = obj.iterMax;
                    fprintf("\n\nMax iteration count reached\n")
                    break
                elseif mod(iter, 10) == 0
                    if obj.rampB(1) == 1
                        if obj.beta < obj.rampB(2)
                            if obj.beta < 1
                                obj.beta = obj.beta + 0.3;
                            else
                                obj.beta = obj.beta * obj.rampB(3);
                            end
                            if obj.beta >= obj.rampB(2)
                                obj.logs.('Heavi_done') = iter;
                            end
                        end
                    end
                    if obj.p < obj.rampPQ(1)
                        obj.p = obj.p + obj.rampPQ(2);
                        obj.q = obj.q + obj.rampPQ(2);
                        if obj.p >= obj.rampPQ(1)
                            obj.logs.('SIMP_done') = iter;
                            simp_iter = iter;
                        end
                    end
                end
                obj = init(obj, x);
                obj = newt(obj);
                [obj.g0(iter), dg0, obj.gc(iter, :), dgc, obj.cp] = funcEval(obj, x);
                if iter == 1
                    s = abs(100/obj.g0(1));
                end
                [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(obj.ncon, obj.nel, iter, x, zeros(obj.nel, 1), ones(obj.nel, 1), ...
                                                            xold1, xold2, s*obj.g0(iter), s*dg0, obj.gc(iter, :)', dgc, low, upp, a0, a1, c, d, obj.mmaVals);
                xold2 = xold1;
                xold1 = x;
  
                obj.dx(iter) = norm(xmma - x, inf);
                if iter > 1000
                    fconv = norm(diff(obj.g0(iter-10:iter)), inf);
                else
                    fconv = inf;
                end
                criteria = (fconv > obj.ftol && obj.dx(iter) > obj.xtol) || any(obj.gc(iter, :) > 0) || obj.p < obj.rampPQ(1) || obj.rampB(1) == 2;

                if obj.rampB(1) == 2 && (obj.dx(iter) < obj.xtol*10 || iter > simp_iter + 50)
                    obj.rampB(1) = 1;
                    obj.logs.('Heavi_on') = iter;
                end

                x = xmma;
                rho = he(obj, (obj.Z*x + obj.Zp));
                if obj.plots
                    if iter == 1
                        obj.design = figure;
                    end
                    plotFigs(obj, rho, 1);
                end
                fprintf("Opt iter: %i\n", iter)
                if obj.stressCon
                    fprintf("  g0: %.2g, g1: %.2g, g2: %.2g, dx: %.2g\n", [s*obj.g0(iter), obj.gc(iter, 1), obj.gc(iter, 2), obj.dx(iter)])
                else
                    fprintf("  g0: %.2g, g1: %.2g, dx: %.2g\n", [s*obj.g0(iter), obj.gc(iter, 1), obj.dx(iter)])
                end
            end
            obj.dx = obj.dx(1:iter);
            obj.g0 = s*obj.g0(1:iter);
            obj.gc = obj.gc(1:iter, :);
            if obj.plots
                plotFigs(obj, rho, 0);
            end
        end

        function [g0, dg0, gc, dgc, cp] = funcEval(obj, x)
            dgt0dx = zeros(1, obj.nel);
            dR1dx = zeros(obj.endof*obj.nel, 1);  
            dR2dx = zeros(obj.tgp, 1);  
            dgt0dep = zeros(obj.tgp, 1);
            dR1dep = zeros(obj.endof*obj.tgp, 1);
            ap = sparse(obj.disp(:, 1), 1, obj.a(obj.disp(:, 1)), obj.ndof, 1);
            
            sigb = zeros(obj.nel, 1); 
            dsigbdx = zeros(obj.nel, 1);
            dsigbdep = zeros(obj.tgp, 1);
            dsigbda = zeros(obj.nel, obj.endof);
            
            [x, dxH] = he(obj, (obj.Z*x + obj.Zp));
            dgam = obj.p*(1-obj.del)*x.^(obj.p-1);
            dphi = obj.q*(1-obj.dels)*x.^(obj.q-1);
            th = (dgam.*obj.phi - dphi.*obj.gam)./obj.phi.^2;
            for el = 1:obj.nel
                Kte = zeros(obj.endof);
                eix = obj.edof(el, :);

                for gp = 1:obj.ngp
                    B = obj.Bgp(3*(gp-1)+1:3*gp, :); J = obj.detJ(gp); % [B, J] = NablaB(obj, gp, el);
                    ix = obj.ngp*(el-1) + gp;
                    ixM =  4*obj.ngp*(el-1) + (gp-1)*4 + 1:4*obj.ngp*(el-1) + gp*4;
                    k0 = obj.ep(ix)*obj.sigy0^2/(obj.sigy0 + obj.H*obj.ep(ix) + obj.Kinf*(1-exp(-obj.xi*obj.ep(ix))));
                    dDsdx = (dgam(el)*obj.Ds(ixM, :) - th(el)*k0*obj.Ds(ixM, :)*obj.P*obj.Ds(ixM, :))/obj.gam(el);
                    Kte = Kte + B'*dDsdx(([1 2 4]), [1 2 4])*B*J*obj.t;

                    depstdx = obj.eps(ix, :)*(dDsdx*obj.P*obj.Ds(ixM, :)...
                              + obj.Ds(ixM, :)*obj.P*dDsdx)*obj.eps(ix, :)';
                    dR2dx(ix) = dphi(el)*(obj.sigy0 + obj.H*obj.ep(ix) + obj.Kinf*(1-exp(-obj.xi*obj.ep(ix))))...
                                - obj.sigy0/2/sqrt(obj.phi(el)^2*obj.epst(ix))*depstdx;

                    Kh = B'*obj.dDsdep(ixM([1 2 4]),[1 2 4])*B*J*obj.t;
                    dR1depe = Kh*obj.a(eix);
                    dgt0dep(ix) = -ap(eix)'*dR1depe;
                    dR1dep(obj.endof*(ix-1)+1:obj.endof*ix) = dR1depe; 
                    
                    if obj.stressCon
                        sigP = obj.sig(ix, :)*obj.P;
                        sigb(el) = sigb(el) + sigP*obj.sig(ix, :)'/obj.ngp;
                        dsigbdx(el) = dsigbdx(el) + obj.sigy0^2*sigP*dDsdx*obj.eps(ix, :)'/obj.ngp;
                        dsigbda(el, :) = dsigbda(el, :) + obj.sigy0^2*sigP*obj.Dt(ixM, [1 2 4])*B/obj.ngp;
                        dsigbdep(ix) = obj.sigy0^2*trace(obj.P*obj.Ds(ixM, :)*obj.eps(ix, :)'*obj.eps(ix, :)*obj.dDsdep(ixM, :))/obj.ngp;
                    end
                end
                dR1dxe = Kte*obj.a(eix);
                dgt0dx(el) = -ap(eix)'*dR1dxe;
                dR1dx(obj.endof*(el-1)+1:obj.endof*el) = dR1dxe; 
                
                if obj.stressCon
                    sigb(el) = obj.sigy0*sqrt(sigb(el));
                    dsigbdep(obj.ngp*(el-1)+1:obj.ngp*el) = dsigbdep(obj.ngp*(el-1)+1:obj.ngp*el)/sigb(el)/obj.phi(el)/obj.sigm(el);
                end
            end
            dgt0da = -obj.a(obj.disp(:, 1))'*obj.K(obj.disp(:, 1), obj.fdof);

            dR1dx = sparse(reshape(obj.edof', [], 1), repelem((1:obj.nel)', obj.endof), dR1dx, obj.ndof, obj.nel);
            dR2dx = sparse((1:obj.tgp)', repelem((1:obj.nel)', obj.ngp), dR2dx, obj.tgp, obj.nel);
            dR1dep = sparse(reshape(repelem(obj.edof', 1, obj.ngp), [], 1), repelem((1:obj.tgp)', obj.endof), dR1dep, obj.ndof, obj.tgp);

            pgp = find(obj.ep);
            lamt = -dgt0da/obj.K(obj.fdof, obj.fdof);
            idR2dep = diag(1./obj.dR2dep(pgp));
            mut = -dgt0dep(pgp)'*idR2dep - lamt*dR1dep(obj.fdof, pgp)*idR2dep;

            g0 = -obj.a(obj.disp(:, 1))'*obj.R1(obj.disp(:, 1));
            dg0 = dxH'.*obj.Z'*(dgt0dx + lamt*dR1dx(obj.fdof, :) + mut*dR2dx(pgp, :))';
            dg0(obj.fixDens) = 0;

            g1 = x'*obj.A/obj.Amax - 1;
            dg1 = (dxH'.*obj.Z'*obj.A/obj.Amax)';
            dg1(obj.fixDens) = 0;

            if ~obj.stressCon
                gc = g1;
                dgc = dg1;
                cp = [];
            else
                sigp = norm(sigb./obj.phi./obj.sigm, obj.pnm);
                dgt2dx = obj.cp*(sigb./obj.phi./obj.sigm/sigp).^(obj.pnm - 1).*(dsigbdx./sigb.*obj.phi - sigb.*dphi)./obj.phi.^2./obj.sigm;
                dgt2da = obj.cp*(sigb./obj.phi./obj.sigm/sigp).^(obj.pnm - 1).*dsigbda./sigb./obj.phi./obj.sigm;
                dgt2da = accumarray(reshape(obj.edof', [], 1), reshape(dgt2da', [], 1), [obj.ndof, 1]);
                dgt2dep = obj.cp*(repelem(sigb./obj.phi./obj.sigm, obj.ngp)/sigp).^(obj.pnm - 1).*dsigbdep;

                nut = -dgt2da(obj.fdof)'/obj.K(obj.fdof, obj.fdof);
                xit = -dgt2dep(pgp)'*idR2dep - nut*dR1dep(obj.fdof, pgp)*idR2dep;

                g2 = obj.cp*sigp - 1;
                dg2 = (dxH'.*obj.Z'*(dgt2dx' + nut*dR1dx(obj.fdof, :) + xit*dR2dx(pgp, :))')';
                dg2(obj.zeroGrad) = 0;
                gc = [g1; g2];
                dgc = [dg1; dg2];

                cp = obj.ca*max(sigb./obj.phi./obj.sigm)/sigp + (1 - obj.ca)*obj.cp;
            end
        end

        function [Z, Zp] = filterMatrix(obj, le, re, pc)
            if obj.filtOn
                ec = [obj.ex(:, 1) + le/2, obj.ey(:, 1) + le/2];
                I = zeros(obj.nel*(2*re)^2, 3);
                [x, y] = meshgrid(-re:re, -re:re);
                weights = max(0, 1 - sqrt(x.^2 + y.^2)/re);
                sw = sum(weights(:));
                r0 = le*re;
                Zp = zeros(obj.nel, 1);
                i = 0;
                for ii = 1:obj.nel
                    r = vecnorm(ec - ec(ii, :), 2, 2);
                    ix = find(r - r0 < 0);
                    ixn = length(ix);
                    w = 1-r(ix)/r0;
                    I(i+1:i+ixn, :) = [ii*ones(ixn, 1), ix, w/sw];
                    i = i + ixn;
                    if ~isnan(pc)
                        rp = vecnorm(pc(:, 1:2) - ec(ii, :), 2, 2);
                        ixp = find(rp - r0 < 0);
                        Zp(ii) = pc(ixp, 3)'*(1-rp(ixp)/r0)/sw;
                    end
                end
                Z = sparse(I(1:i, 1), I(1:i, 2), I(1:i, 3), obj.nel, obj.nel);
            else
                Z = speye(obj.nel);
                Zp = zeros(obj.nel, 1);
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
                    obj.sig1N(:,1:4) = obj.sig;
                end
                n = n + 1;
            end
            obj.sig1N(:,5:8) = obj.sig;
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
            if nargin == 3
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
                depstdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
                dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
                drdep = (obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(depstdDs*dDsdep));
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
            depstdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
            dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
            drdep =  phi*(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(depstdDs*dDsdep));
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
                x = he(obj, (obj.Z*x + obj.Zp));
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
                if isempty(obj.design)
                    obj.design = figure;
                end
                figure(obj.design);
                clf;
                colormap(flipud(gray(256)));
                patch(obj.ex', obj.ey', x, ...
                    'EdgeColor', 'none');
                axis equal off;
                axis(obj.axi);
                drawnow;
            else
                cosT = zeros(obj.nel,1);
                Hs = zeros(obj.nel, 1);
                for el = 1:obj.nel
                    ix = obj.ngp*(el-1)+1:el*obj.ngp;
                    cosT(el) = trace(obj.sig1N(ix,1:4)*obj.sig1N(ix,5:8)')/(vecnorm(obj.sig1N(ix,1:4),2,2)'*vecnorm(obj.sig1N(ix,5:8),2,2));
                    Hs(el) = sqrt(obj.sigy0^2*trace(obj.sig(ix, :)*obj.P*obj.sig(ix, :)')/obj.ngp);
                end
                cosT(x<0.5) = 0;
                Hc = Hs./(obj.phi.*obj.sigy(1:4:end));
                Hc(x<0.5) = 0;

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

                % norm(maxk(Hs./obj.phi./obj.sigm, 10), 8)
                % sigp = norm(Hs./obj.phi./obj.sigm, obj.pnm+8)
                % Hs(~ismember(1:numel(Hs), idx)) = 0;
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
                legend('Location', 'south', 'FontSize', 12);
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
                    ylim([round(min(-obj.g0(10:end))-6, -1) round(max(-obj.g0(10:end)+5), -1)]);
                    ax = gca;
                    ax.YColor = ax.XColor;

                    yyaxis right;
                    hold on;
                    plot((10:iter)', obj.gc(10:end, 1), 'k', 'LineWidth', 2);
                    if obj.stressCon
                        plot((10:iter)', obj.gc(10:end, 2), 'b', 'LineWidth', 2);
                        ylabel('g_{1,2} constraint');
                        legend('Stiffness', 'Volume', 'Stress');
                    else
                        ylabel('g1 constraint');
                        legend('Stiffness', 'Volume');
                    end
                    hold off;
                    ylim([-0.05 0.05]);
                    ax = gca;
                    ax.YColor = ax.XColor;
                    xlabel('Iteration');
                    title('Convergence Plot');
                    grid on;

                    figure;
                    plot((10:iter)', obj.dx(10:end), 'b', 'LineWidth', 2);
                    ylim([0 ceil(max(obj.dx)*11)/11]);
                    xlabel('Iteration');
                    ylabel('dx');
                    title('Design Change Plot');
                    grid on;
                end
            end
        end

        function plotGrads(obj, dg0, dgc)
            figure;
            patch(obj.ex', obj.ey', dg0, ...
                    'EdgeColor', 'none');
            colormap("hot");
            colorbar;
            title('dg0');
            axis equal off;
            axis(obj.axi);

            figure;
            patch(obj.ex', obj.ey', dgc(2, :), ...
                    'EdgeColor', 'none');
            colormap("hot");
            colorbar;
            title('dgc');
            axis equal off;
            axis(obj.axi);
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
            fields = {'ex', 'ey', 'ngp', 'nel', 'sig', 'P', 'sigy0', 'sig1N', 'ep', 'g0', 'gc', 'dx', 'beta', 'logs', 'sigy'};
            for i = 1:numel(fields)
                if isfield(in, fields{i}) || isprop(in, fields{i})
                    out.(fields{i}) = in.(fields{i});
                elseif strcmp(fields{i}, 'dx')
                    out.(fields{i}) = ones(length(in.g0), 1);
                elseif strcmp(fields{i}, 'beta')
                    out.(fields{i}) = 10;
                elseif strcmp(fields{i}, 'logs')
                    out.(fields{i}) = [];
                elseif strcmp(fields{i}, 'sigy')
                    out.(fields{i}) = ones(in.nel*in.ngp, 1)*in.sigy0;
                end
            end
        end

        function bc = addBC(bc, ly, le, ndof)
            nR = ly/le + 1;
            fix = [ndof/nR*(1:nR)'; ndof/nR*(1:nR)'-1];
            bc = [bc; [fix zeros(2*nR,1)]];
        end

        function pc = padding(lx, ly, le, wx, wy, re, loadcase, pad)
            if loadcase == 3 && pad
                XC = {le/2:le:lx+le*(re-0.5), lx+le/2:le:lx+le*(re-0.5), lx-le/2:-le:wy+le/2,...
                     wy+le/2:le:wy+le*(re-0.5), wy-le/2:-le:-le*(re-0.5), -le/2:-le:-le*(re-0.5)};
                YC = {-le/2:-le:-le*(re-0.5), le/2:le:wx+le*(re-0.5), wx+le/2:le:wx+le*(re-0.5),...
                     wx+le*(re+0.5):le:ly+le*(re-0.5), ly+le/2:le:ly+le*(re-0.5), ly-le/2:-le:-le*(re-0.5)};
                lenX = cellfun(@numel, XC);
                lenY = cellfun(@numel, YC);
                pc = zeros(lenX*lenY', 3);
                nP = 1;
                for i = 1:6
                    nX = numel(XC{i});
                    nY = numel(YC{i});
                    pc(nP:nP+nX*nY-1, 1) = repelem(XC{i}, nY);
                    pc(nP:nP+nX*nY-1, 2) = repmat(YC{i}, 1, nX);
                    nP = nP + numel(XC{i})*numel(YC{i});
                end
                if 3e-3/le < 2
                    bl = 2*le;
                else
                   bl = 3e-3;
                end
                pc(:, 3) = double(pc(:,1) > lx  & abs(pc(:, 2) - wx/2) < bl + 2*le);
            else
                pc = NaN;
            end
        end
        
        function drawDesign(sol, val, x, plots)
            sol = sol.assignVar(val, sol);
            rho = sol.he(sol.Z*x + sol.Zp);
            sol.phi = sol.dels + (1-sol.dels)*rho.^sol.q;
            plotFigs(sol, rho, 1);
            if ~plots
                plotFigs(sol, rho, 0);
            end
        end
        
        function drawMultipleDesigns(JOB)
            path = fullfile("batch", sprintf("JOB_%s", JOB));
            jobs = dir(fullfile(path, "**", "*.mat"));
            jobs = jobs(~strcmp({jobs.name}, "input.mat"));
            for i=1:length(jobs)
                load(fullfile(jobs(i).folder, jobs(i).name), "params", "val", "x");
                if i == 1
                    sol = Solver(params);
                end
                sol.drawDesign(sol, val, x, 1);
                title(jobs(i).name, 'Interpreter', 'none')
            end
        end
    end
end