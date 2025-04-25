clc, clear, close all

% Input parameters
params.le = 0.05;
params.lx = 1;
params.ly = 1;
params.Vf = 0.25;

params.t = 1;
params.ngp = 4;

params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy0; %extra terms linHard (sat. stress)
params.xi = 1e-3; %extra terms linHard (sat. exp)
params.rtol = 1e-5;

E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;
params.v21 = v; params.v32 = v; params.v31 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;

params.DP = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

params.N = 5;
params.R1tol = 1e-4;
params.disp = [2 -4e-4;
               4 -4e-4;
               6 -4e-4]; % displacement [nodes total-size]

% Optimization Parameters
params.re = 2; % Elements in radius
params.p = 3;
params.q = 2;
params.del = 1e-9;
params.ncon = 1; % Nr of constraints
params.xTol = 1e-4;
params.iterMax = 750;

sol = Solver(params);

%% Opt test
x = ones(sol.nel, 1);
sol = optimizer(sol, x);

%% Mesh
patch(sol.ex', sol.ey', 1)

%% Newton-Raphson
sol.gam = ones(sol.nel, 1);
sol.phi = ones(sol.nel, 1);
sol = newt(sol);
plotFigs(sol, ones(sol.nel, 1), 1)

figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);
dof = sol.ndof-3;
fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

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

%% Finite diff
sol1 = Solver(params);
sol2 = Solver(params);
h = 1e-6;

% c = [0.3 0.5 0.2 0.7 0.9]';
% x = repmat(c, sol.nel/5, 1);
x = 0.5*ones(sol.nel, 1);
x1 = x;
x2 = x;
el = 11;

x1(el) = x1(el) - h;
x2(el) = x2(el) + h;
[~, ~, dg, ~] = optimizer(sol, x);
[~, g1, ~, ~] = optimizer(sol1, x1);
[~, g2, ~, ~] = optimizer(sol2, x2);

dgf = (g2-g1)/2/h;
fprintf("\nDiff: %.5g \ndg: %.5g \ndgf: %.5g", [dgf-dg(el), dg(el), dgf])

