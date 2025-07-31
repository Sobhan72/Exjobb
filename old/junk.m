%% Mesh
patch(sol.ex', sol.ey', ones(sol.nel, 1));
colorbar;
colormap jet;

%% Model validation
params.DP = 0;
sol = Solver(params);
sol = newt(sol);
dof = sol.ndof/2+1;
disp1 = sol.a(dof);
sol.plotFigs

params.DP = 1;
params.N = 5;
sol = Solver(params);
sol = newt(sol);
disp2 = sol.a(dof);
sol.plotFigs
fprintf("\nDisp1: %.5g \nDisp2: %.5g \nDiff: %.5g%% \n", [disp1, disp2, (1-disp2/disp1)*100])

%% Hill material model test
N = 5; 
reverse = 0;

epse = zeros(N,1); % eps eff
sige = zeros(N,1); % sig eff
siggp = [0;0;0;0];
epsgp = [0;0;0;0];
deps = [1e-4, 0, 0, 1e-4]'*10;
sigegp = 0;
sigygp = 360e6;
Dsgp = sol.De;
epgp = 0;

R = 2/3*[2/3, -1/3, -1/3,  0;
        -1/3,  2/3, -1/3,  0;
        -1/3, -1/3,  2/3,  0;
         0,     0,    0,   2];

for n = 1:N
    if n == 30 && reverse
        deps = -deps;
    elseif n == 40 && reverse
        deps = -deps;
    end
    fprintf("Load step: %i \n", n)
    epsgp = epsgp + deps;
    epse(n) = sqrt(epsgp'*R*epsgp);
    depsm = deps + [0;0;0;1]*deps(4);
    epsgpm = epsgp + [0;0;0;1]*epsgp(4);
    sigtr = sol.De*depsm + siggp;
    if sqrt(sol.sigy0^2*sigtr'*sol.P*sigtr) > sigygp
        if sol.DP
            [Dtgp, siggp, Dsgp, epgp] = DPMat(sol, epsgpm, Dsgp, epgp);
        else
            [Dtgp, siggp, sigygp] = EPMat(sol, sigtr, sigygp);
        end
    else
        siggp = sigtr;
        Dtgp = sol.De;
    end
    sige(n) = sqrt(sol.sigy0^2*siggp'*sol.P*siggp);
end

figure;
plot([0; epse], [0; sige]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
title("Hill Deformation Model")
grid on;

%% FEM model test
epse = zeros(sol.N, 1);
sige = zeros(sol.N, 1);
gp = 12;
R = 2/3*[2/3, -1/3, -1/3,  0;
        -1/3,  2/3, -1/3,  0;
        -1/3, -1/3,  2/3,  0;
         0,     0,    0,   1/2];

for n = 1:sol.N
    fprintf("Load step: %i \n", n);
    bcD = sol.disp;
    Nr = 0;
    while norm(sol.r1) > sol.r1tol || Nr == 0
        Nr = Nr + 1;
        sol = FEM(sol, bcD);
        bcD(:, 2) = bcD(:, 2)*0;
        fprintf("  Nr: %i, r1: %4.2g \n", [Nr, norm(sol.r1)]);
    end
    sol.eps = sol.epsi; sol.sig = sol.sigi; sol.ep = sol.epi; sol.Ds = sol.Dsi; sol.sigy = sol.sigyi;
    epse(n) = sqrt(sol.eps(gp, :)*R*sol.eps(gp, :)');
    sige(n) = sqrt(sol.sigy0^2*sol.sig(gp, :)*sol.P*sol.sig(gp, :)');
end
sol.plotFigs
figure;
plot([0; epse], [0; sige]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
title("Effective stress-strain")
grid on;


%% Finite diff dx
h = 1e-6;
sol1 = Solver(params);
sol2 = Solver(params);

% c = [0.3 0.5 0.2 0.7 0.9]';
% x = repmat(c, sol.nel/5, 1);
x = ones(sol.nel, 1);
x1 = x;
x2 = x;
elc = 72;
x1(elc) = x1(elc) - h;
x2(elc) = x2(elc) + h;

sol = initOpt(sol, x);
sol = newt(sol);
[~, ~, ~, ~, dg] = funcEval(sol, x);

sol1.ep = sol.ep;
sol1.eps = sol.eps;
sol1.a = sol.a;
sol2.ep = sol.ep;
sol2.eps = sol.eps;
sol2.a = sol.a;

sol1.gam = sol1.del + (1-sol1.del)*x1.^sol1.p;
sol1.phi = sol1.del + (1-sol1.del)*x1.^sol1.q;

sol2.gam = sol2.del + (1-sol2.del)*x2.^sol2.p;
sol2.phi = sol2.del + (1-sol2.del)*x2.^sol2.q;


tripf1 = zeros(sol.nel*sol.endof, 2);
tripf2 = zeros(sol.nel*sol.endof, 2);
ii = 0;
for el = 1:sol.nel
    gam1 = sol1.gam(el);
    gam2 = sol2.gam(el);
    phi1 = sol1.phi(el);
    phi2 = sol2.phi(el);

    fein1 = zeros(sol.endof, 1);
    fein2 = zeros(sol.endof, 1);

    for gp = 1:sol.ngp
        ix = sol.ngp*(el-1) + gp;
        if ix == 285
            1;
        end
        ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        if sol1.ep(ix) ~= 0
            sol1.Ds(ixM, :) = gam1*sol.X*diag(1./diag(eye(4) + gam1/phi1*sol.sigy0^2/(sol.sigy0+sol.H*sol1.ep(ix))*sol1.ep(ix)*sol.Gam))*sol.X';
            
            epst = 1/phi1^2*(sol1.eps(ix, :)*sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*sol1.eps(ix, :)');
            sige = sol.sigy0 + sol.H*sol1.ep(ix);
            r1 = phi1*(sige - sol.sigy0*sqrt(epst));
            % detdDs = 1/phi1^2*(2*sol.P*sol1.Ds(ixM, :)*(sol1.eps(ix, :)'*sol1.eps(ix, :)));
            % dDsdep = 1/phi1*(-sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*(sol.sigy0^2*(sol.sigy0+sol.H*sol1.ep(ix)-sol.H*sol1.ep(ix))/(sol.sigy0+sol.H*sol1.ep(ix))^2));
            % drdep =  phi1*(sol.H - sol.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            % drdeps = 1/phi1*(-sol.sigy0/sqrt(epst)*sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*sol1.eps(ix, :)');
            % depdeps = -drdeps/drdep;
            % sol1.Dt(ixM, :) = sol1.Ds(ixM, :) + dDsdep*sol1.eps(ix, :)'*depdeps';
        else
            % sol1.Dt(ixM, :) = gam1*sol.De;
            sol1.Ds(ixM, :) = gam1*sol.De;
        end

        if sol2.ep(ix) ~= 0
            sol2.Ds(ixM, :) = gam2*sol.X*diag(1./diag(eye(4) + gam2/phi2*sol.sigy0^2/(sol.sigy0+sol.H*sol2.ep(ix))*sol2.ep(ix)*sol.Gam))*sol.X';

            epst = 1/phi2^2*(sol2.eps(ix, :)*sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*sol2.eps(ix, :)');
            sige = sol.sigy0 + sol.H*sol2.ep(ix);
            r2 = phi2*(sige - sol.sigy0*sqrt(epst));
            % detdDs = 1/phi2^2*(2*sol.P*sol2.Ds(ixM, :)*(sol2.eps(ix, :)'*sol2.eps(ix, :)));
            % dDsdep = 1/phi2*(-sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*(sol.sigy0^2*(sol.sigy0+sol.H*sol2.ep(ix)-sol.H*sol2.ep(ix))/(sol.sigy0+sol.H*sol2.ep(ix))^2));
            % drdep =  phi2*(sol.H - sol.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            % drdeps = 1/phi2*(-sol.sigy0/sqrt(epst)*sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*sol2.eps(ix, :)');
            % depdeps = -drdeps/drdep;
            % sol2.Dt(ixM, :) = sol2.Ds(ixM, :) + dDsdep*sol2.eps(ix, :)'*depdeps';
        else
            % sol2.Dt(ixM, :) = gam2*sol.De;
            sol2.Ds(ixM, :) = gam2*sol.De;
        end

        [B, J] = NablaB(sol, gp, el);
        sig1 = sol1.Ds(ixM, :)*sol1.eps(ix, :)';
        dsig = norm(sol.sig(ix, :) - sig1');
        % if dsig>1
        %     ix
        % end
        sig2 = sol2.Ds(ixM, :)*sol2.eps(ix, :)';
        fein1 = fein1 + B'*sig1([1 2 4])*J*sol.t;
        fein2 = fein2 + B'*sig2([1 2 4])*J*sol.t;
        dR2dx1(ix) = r1;
        dR2dx2(ix) = r2;
    end
    tripf1((el-1)*sol.endof+1:el*sol.endof, :) = [sol.edof(el, :)', fein1];
    tripf2((el-1)*sol.endof+1:el*sol.endof, :) = [sol.edof(el, :)', fein2];
end

fin1 = sparse(tripf1(:, 1), 1, tripf1(:, 2), sol.ndof, 1);
fin2 = sparse(tripf2(:, 1), 1, tripf2(:, 2), sol.ndof, 1);
sol1.R1 = fin1;
sol2.R1 = fin2;
sol1.K = sol.K;
sol2.K = sol.K;

% [gf1, ~, ~, ~]  = funcEval(sol1, x1);
% [gf2, ~, ~, ~]  = funcEval(sol2, x2);

dgf = (dR2dx2 - dR2dx1)/(2*h);
% dgf = (gf2-gf1)/2/h;
% dff = (sol2.R1 - sol1.R1)/(2*h);
% dDs = (sol2.Ds(4*sol.ngp*(elc-1) + 1:4*sol.ngp*(elc-1) + 16, :) - sol1.Ds(4*sol.ngp*(elc-1) + 1:4*sol.ngp*(elc-1) + 16, :))/(2*h);
% fprintf("\nDiff: %.5g \ndg: %.5g \ndgf: %.5g", [dgf-dg(elc), dg(elc), dgf])

%% Finite diff dep
h = 1e-8;
sol1 = Solver(params);
sol2 = Solver(params);

x = 0.8*ones(sol.nel, 1);
[sol, ~, ~, ~, ~, ~, dg] = optimizer(sol, x);

ep1 = sol.ep;
ep2 = sol.ep;
gpc = 17;
ep1(gpc) = ep1(gpc) - h;
ep2(gpc) = ep2(gpc) + h;

sol1.eps = sol.eps;
sol2.eps = sol.eps;

tripf1 = zeros(sol.nel*sol.endof, 2);
tripf2 = zeros(sol.nel*sol.endof, 2);
for el = 1:sol.nel
    gam1 = sol.gam(el);
    gam2 = sol.gam(el);
    phi1 = sol.phi(el);
    phi2 = sol.phi(el);

    fein1 = zeros(sol.endof, 1);
    fein2 = zeros(sol.endof, 1);

    for gp = 1:sol.ngp
        ix = sol.ngp*(el-1) + gp;
        ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        if sol.ep(ix) ~= 0
            sol1.Ds(ixM, :) = gam1*sol.X*diag(1./diag(eye(4) + gam1/phi1*sol.sigy0^2/(sol.sigy0+sol.H*ep1(ix))*ep1(ix)*sol.Gam))*sol.X';
            
            epst = 1/phi1^2*(sol1.eps(ix, :)*sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*sol1.eps(ix, :)');
            detdDs = 1/phi1^2*(2*sol.P*sol1.Ds(ixM, :)*(sol1.eps(ix, :)'*sol1.eps(ix, :)));
            dDsdep = 1/phi1*(-sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*(sol.sigy0^2*(sol.sigy0+sol.H*ep1(ix)-sol.H*ep1(ix))/(sol.sigy0+sol.H*ep1(ix))^2));
            drdep =  phi1*(sol.H - sol.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            drdeps = 1/phi1*(-sol.sigy0/sqrt(epst)*sol1.Ds(ixM, :)*sol.P*sol1.Ds(ixM, :)*sol1.eps(ix, :)');
            depdeps = -drdeps/drdep;
            sol1.Dt(ixM, :) = sol1.Ds(ixM, :) + dDsdep*sol1.eps(ix, :)'*depdeps';
        else
            sol1.Dt(ixM, :) = gam1*sol.De;
            sol1.Ds(ixM, :) = gam1*sol.De;
        end
        

        if sol.ep(ix) ~= 0
            sol2.Ds(ixM, :) = gam2*sol.X*diag(1./diag(eye(4) + gam2/phi2*sol.sigy0^2/(sol.sigy0+sol.H*ep2(ix))*ep2(ix)*sol.Gam))*sol.X';
    
            epst = 1/phi2^2*(sol2.eps(ix, :)*sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*sol2.eps(ix, :)');
            detdDs = 1/phi2^2*(2*sol.P*sol2.Ds(ixM, :)*(sol2.eps(ix, :)'*sol2.eps(ix, :)));
            dDsdep = 1/phi2*(-sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*(sol.sigy0^2*(sol.sigy0+sol.H*ep2(ix)-sol.H*ep2(ix))/(sol.sigy0+sol.H*ep2(ix))^2));
            drdep =  phi2*(sol.H - sol.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
            drdeps = 1/phi2*(-sol.sigy0/sqrt(epst)*sol2.Ds(ixM, :)*sol.P*sol2.Ds(ixM, :)*sol2.eps(ix, :)');
            depdeps = -drdeps/drdep;
            sol2.Dt(ixM, :) = sol2.Ds(ixM, :) + dDsdep*sol2.eps(ix, :)'*depdeps';
        else
            sol2.Dt(ixM, :) = gam2*sol.De;
            sol2.Ds(ixM, :) = gam2*sol.De;  
        end

        [B, J] = NablaB(sol, gp, el);
        sig1 = sol1.Ds(ixM, :)*sol1.eps(ix, :)';
        sig2 = sol2.Ds(ixM, :)*sol2.eps(ix, :)';
        fein1 = fein1 + B'*sig1([1 2 4])*J*sol.t;
        fein2 = fein2 + B'*sig2([1 2 4])*J*sol.t;
    end
    tripf1((el-1)*sol.endof+1:el*sol.endof, :) = [sol.edof(el, :)', fein1];
    tripf2((el-1)*sol.endof+1:el*sol.endof, :) = [sol.edof(el, :)', fein2];
end

fin1 = sparse(tripf1(:, 1), 1, tripf1(:, 2), sol.ndof, 1);
fin2 = sparse(tripf2(:, 1), 1, tripf2(:, 2), sol.ndof, 1);
sol1.R1 = fin1;
sol2.R1 = fin2;

dgf = (sol2.R1 - sol1.R1)/(2*h);



%%


for n = 1:sol.N
    fprintf("Load step: %i \n", n);
    bcD = sol.disp;
    Nr = 0;
    while norm(sol.r1) > sol.r1tol || Nr == 0
        Nr = Nr + 1;
        sol = FEM(sol, bcD);
        bcD(:, 2) = bcD(:, 2)*0;
        fprintf("  Nr: %i, r1: %4.2g \n", [Nr, norm(sol.r1)]);
    end
    sol.eps = sol.epsi; sol.sig = sol.sigi; sol.ep = sol.epi; sol.Ds = sol.Dsi; sol.sigy = sol.sigyi;
end

gp = 12;
p = 6;
q = 5;
d = 1e-9;
x = 0.5;
h = 1e-8;

k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';

dgam = p*(1-d)*x^(p-1); 
dphi = q*(1-d)*x^(q-1);
th = (dgam*phi-dphi*gam)/(phi)^2;
dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;
dPdx = -2/phi^3*sol.P;
depstdx = sol.eps(gp, :)*dDsdx*sol.P/phi^2*Ds*sol.eps(gp, :)' - sol.eps(gp, :)*Ds*sol.P*2*dphi/phi^3*Ds*sol.eps(gp, :)' + sol.eps(gp, :)*Ds*sol.P/phi^2*dDsdx*sol.eps(gp, :)';
dR2dx = dphi*(sol.sigy0+sol.H*sol.ep(gp))-dphi*sol.sigy0*sqrt(epst)-phi*sol.sigy0/2/sqrt(epst)*depstdx;

x = 0.5 - h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';
r1 = phi*(sol.sigy0 + sol.H*sol.ep(gp)) - phi*sol.sigy0*sqrt(epst);

x = 0.5 + h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';
r2 = phi*(sol.sigy0 + sol.H*sol.ep(gp)) - phi*sol.sigy0*sqrt(epst);

dr = (r2-r1)/(2*h);

dr - dR2dx




k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
dgam = p*(1-d)*x^(p-1); 
dphi = q*(1-d)*x^(q-1);
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
h = 1e-8;

th = (dgam*phi-dphi*gam)/(phi)^2;

dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;

x=0.5-h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;

x=0.5+h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds2 = gam*sol.De*V*sol.De;

dDs = (Ds2-Ds)/(2*h);
dDsdx-dDs;


tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        dgam = p*(1-d)*x^(p-1); 
        dphi = q*(1-d)*x^(q-1);
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        
        th = (dgam*phi-dphi*gam)/(phi)^2;
        
        dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;

    ke = ke + B'*dDsdx([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
Kt = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g0t = params.disp(:, 2)'*Kt(params.disp(:, 1),:)*disp2;

x = 0.5 - h;
tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        Ds = gam*sol.De*V*sol.De;

    ke = ke + B'*Ds([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
K1 = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g01 = params.disp(:, 2)'*K1(params.disp(:, 1),:)*disp2;

x = 0.5 + h;
tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        Ds = gam*sol.De*V*sol.De;

    ke = ke + B'*Ds([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
K2 = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g02 = params.disp(:, 2)'*K2(params.disp(:, 1),:)*disp2;

dg = (g02-g01)/(2*h);

g0t-dg;

%% Func.
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

Dti = [];

eps = [eps(1) eps(4)/2 0;
    eps(4)/2 eps(2) 0;
    0 0 eps(3)];

depdeps = [depdeps(1) depdeps(4) 0;
    depdeps(4) depdeps(2) 0;
    0 0 depdeps(3)];

dDsdep = tensor(obj, dDsdep);

for ii = indx'
    Dtii = 0;
    for k = 1:3
        for l = 1:3
            Dtii = Dtii + dDsdep(ii(1), ii(2), k, l)*depdeps(ii(3), ii(4))*eps(k,l);
        end
    end
    Dti = [Dti, Dtii];
end

D = Ds + [Dti(1), Dti(2), Dti(3), Dti(4);
    Dti(2), Dti(5), Dti(6), Dti(7);
    Dti(3), Dti(6), Dti(8), Dti(9);
    Dti(4), Dti(7), Dti(9), Dti(10)];
end

function t = tensor(obj, matrix)
t = zeros(3,3,3,3);
t(1,1,1,1) = matrix(1,1);
t(1,1,1,2) = matrix(1,4);
t(1,1,2,1) = matrix(1,4);
t(1,1,2,2) = matrix(1,2);
t(1,1,3,3) = matrix(1,3);

t(1,2,1,1) = matrix(4,1);
t(1,2,1,2) = matrix(4,4);
t(1,2,2,1) = matrix(4,4);
t(1,2,2,2) = matrix(4,2);
t(1,2,3,3) = matrix(4,3);

t(2,1,1,1) = matrix(4,1);
t(2,1,1,2) = matrix(4,4);
t(2,1,2,1) = matrix(4,4);
t(2,1,2,2) = matrix(4,2);
t(2,1,3,3) = matrix(4,3);

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

function Dt = newDt(obj, sige, Ds, eps, delta_ep) %Not correct
sig = Ds*eps;
dfdsig = obj.sig_y0^2/sige*obj.P*sig;
df2dsig2 = obj.sig_y0^2/sige*obj.P*(eye(4)-obj.sig_y0^2/sige^2*sig*sig'*obj.P);
Da = inv(obj.C + delta_ep*df2dsig2);
A = dfdsig'*Da*dfdsig + obj.H;
Dt = Da - 1/A*Da*(dfdsig*dfdsig')*Da;
end

function sig = update_stress(obj, eps)
e = [eps(1:3)-mean(eps(1:3)); 2*eps(4)];
sigkk = 3*obj.Kmod*sum(eps(1:3));
s = obj.G*e;
sig = [s(1:3)+sigkk/3; s(4)];
end

function sig = update_stress_el(obj, eps)
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
sigkk = 3*obj.Kmod*sum(eps(1:3));
s = 2*obj.Ge*e;
sig = [s(1:3)+sigkk/3; s(4)];
end

function stress_eff_out = stress_eff(~, sig)    % Calculate effective strain
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
J2 = 0.5*(s'*s + s(4)*s(4));
stress_eff_out = sqrt(J2*3);
end

function D_ep = Dtan(obj, sig, sig_eff)   %%FEL%%
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
A = obj.H + (obj.sig_y0^2/sig_eff)^2*s'*obj.P*D*obj.P*s; %sig_y = sig_eff
D_p = 1/A*(obj.sig_y0^2/sig_eff)^2*D*obj.P*s*s'*obj.P*D;
D_ep = D-D_p;
end

function [bc, disp] = tempBc(~, ly, le, ndof)
nR = ly/le + 1;
bc = [[ndof/nR*(0:nR-1)'+1; ndof/nR*(0:nR-1)'+2] zeros(2*nR,1)];
disp = [ndof/nR*(1:nR)'-1, ones(nR,1)*3e-3];
end

function Dtf = fdmDt(obj, eps, Ds, ep, gam, phi)
    delta = 1e-9;
    
    Dtf = zeros(4);
    for i = 1:4
        deps = [0; 0; 0; 0];
        deps(i) = delta;
        eps2 = eps - deps;
        eps3 = eps + deps;
        sig2 = DPMat(obj, eps2, Ds, ep, gam, phi);
        sig3 = DPMat(obj, eps3, Ds, ep, gam, phi);
        Dtf(:, i) = (sig3-sig2)/2/delta;
    end
end

function fdtEP(obj, sigtr, sigy, gam, phi)
    Dep = 1e-4;
    U = obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/(sigy + obj.H*Dep)*Dep*obj.Gam))*obj.iX;
    sigt = 1/phi^2*(sigtr'*U'*obj.P*U*sigtr);
    dUdDep = gam/phi*(-U*obj.DeP*U*(obj.sigy0^2*sigy/(sigy + obj.H*Dep)^2));
    dsigtdU = 1/phi^2*(2*obj.P*U*(sigtr*sigtr'));
    drdDep = obj.H - obj.sigy0/(2*sqrt(sigt))*trace(dsigtdU'*dUdDep);
    h = 1e-10;
    Dep1 = Dep + h;
    Dep2 = Dep - h;
    U1 = obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/(sigy + obj.H*Dep1)*Dep1*obj.Gam))*obj.iX;
    U2 = obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/(sigy + obj.H*Dep2)*Dep2*obj.Gam))*obj.iX;
    sigt1 = 1/phi^2*(sigtr'*U1'*obj.P*U1*sigtr);
    sigt2 = 1/phi^2*(sigtr'*U2'*obj.P*U2*sigtr);
    r1 = sigy + obj.H*Dep1 - obj.sigy0*sqrt(sigt1);
    r2 = sigy + obj.H*Dep2 - obj.sigy0*sqrt(sigt2);
    drdDepf = (r1-r2)/2/h;
    fprintf("Diff: %.5g\n", norm(drdDepf-drdDep));
end
function fdtDP(obj, eps, gam, phi)
    ep = 1e-4;
    sige = (obj.sigy0 + obj.H*ep + obj.Kinf*(1-exp(-obj.xi*ep)));
    Ds = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige*ep*obj.Gam))*obj.X';
    detdDs = 1/phi^2*(2*obj.P*Ds*(eps*eps'));
    dDsdep = 1/phi*(-Ds*obj.P*Ds*(obj.sigy0^2*(sige-(obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep))*ep)/sige^2));
    epst = 1/phi^2*(eps'*Ds*obj.P*Ds*eps);
    drdep = (obj.H + obj.Kinf*obj.xi*exp(-obj.xi*ep) - obj.sigy0/(2*sqrt(epst))*trace(detdDs*dDsdep));
    h = 1e-10;
    ep1 = ep + h;
    ep2 = ep - h;
    sige1 = (obj.sigy0 + obj.H*ep1 + obj.Kinf*(1-exp(-obj.xi*ep1)));
    sige2 = (obj.sigy0 + obj.H*ep2 + obj.Kinf*(1-exp(-obj.xi*ep2)));
    Ds1 = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige1*ep1*obj.Gam))*obj.X';
    Ds2 = gam*obj.X*diag(1./diag(eye(4) + gam/phi*obj.sigy0^2/sige2*ep2*obj.Gam))*obj.X';
    epst1 = 1/phi^2*(eps'*Ds1*obj.P*Ds1*eps);
    epst2 = 1/phi^2*(eps'*Ds2*obj.P*Ds2*eps);
    r1 = (sige1 - obj.sigy0*sqrt(epst1));
    r2 = (sige2 - obj.sigy0*sqrt(epst2));
    drdepf = (r1-r2)/2/h;
    fprintf("Diff: %.5g\n", norm(drdepf-drdep));
end
