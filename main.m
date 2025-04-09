clc, clear, close all

% Input parameters
params.le = 0.05;
params.lx = 1;
params.ly = 1;

params.t = 2;
params.ngp = 4;

params.sigy0 = 360e6;
params.H = 10e9;

params.Kinf = 0; %0.3*params.sigy0; %extra terms linHard (sat. stress)
params.del = 1e-3; %extra terms linHard (sat. exp)

params.r2tol = 1e-5;
params.DP = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;
params.v21 = v; params.v32 = v; params.v31 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;

params.N = 10;
params.r1tol = 1e-4;
params.disp = [2 -4e-2;
               4 -4e-2;
               6 -4e-2]; % displacement [nodes total-size]

params.re = 2; % Elements in radius

sol = Solver(params);

%% Mesh
patch(sol.ex', sol.ey', 1)

%% Newton-Raphson
sol = newt(sol);
sol.plotFigs

% figure;
% eldraw2(sol.ex, sol.ey, [1 2 1]);
% hold on
% eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);
% dof = sol.ndof-3;
% fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

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


%% Finite diff
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

p = 6;
q = 5;
d = 1e-9;
x=0.5;

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
g0t = params.disp(:, 2)'*Kt(params.disp(:, 1),:)*ones(sol.ndof, 1)*1e-2;

x = 0.5 - 1e-8;
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
g01 = params.disp(:, 2)'*K1(params.disp(:, 1),:)*ones(sol.ndof, 1)*1e-2;

x = 0.5 + 1e-8;
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
g02 = params.disp(:, 2)'*K2(params.disp(:, 1),:)*ones(sol.ndof, 1)*1e-2;

dg = (g02-g01)/(2e-8);

g0t-dg
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
    % epse(n) = strain_eff(epsgp);
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

function strain_eff_out = strain_eff(eps)    % Calculate effective stress
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
strain_eff_out = sqrt(2/3*(e'*e + e(4)*e(4)));
end
