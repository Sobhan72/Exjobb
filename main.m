clc, clear

% FEM parameters
params.le = 0.005;
params.lx = 0.1;
params.ly = 0.05;
params.Vf = 0.3;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e1;
params.N = 50;
params.disp = -1e-3; % displacement [nodes total-size]

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;
 
params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; %extra terms linHard (sat. stress)
params.xi = 1e3; %extra terms linHard (sat. exp)

params.rtol = 1e-4;
params.DP = 0; % 1 for Deformation plasticity, 0 for Elastoplasticity

% Optimization Parameters
params.re = 3; % Elements in radius
params.filtOn = true;
params.loadcase = 1;
params.p = 1;
params.q = 1;
params.eta = 0.5;
params.beta = 1;
params.rampB = 1; % 0: off, 1: B*1.1, 2: B + 1
params.rampPQ = true;
params.del = 1e-9;
params.dels = 1e-3;
params.ncon = 1; % Nr of constraints
params.xtol = 1e-5;
params.iterMax = 500;

params.print = [1,0,0]; %[Load step, R1, R2] 
params.saveName = "";
sol = Solver(params);

%% Optimization
x = 0.8*ones(sol.nel, 1);
[sol, x] = optimizer(sol, x);
saveData(sol, x, params, "data");

%% Draw Design
clc, clear, close all
load("Anisotrop_x=10_case1_disp=10e-4.mat");
sol = Solver(params);
sol = sol.assignVar(val, sol);
sol.beta = 10; sol.p = 3; sol.q = 2.5; sol.t = 1;
x = sol.He(sol.Z*x);
plotFigs(sol, x, 0);

%% Newton-Raphson
% clear, close all
% load("Anisotrop_x=00_case2_disp=8e-4.mat");
% params.N = 500;
% params.print = [1,1,0];
% params.R1tol = 1e1;
% sol = Solver(params);
% sol = sol.assignVar(val, sol);
% sol.beta = 10; sol.p = 3; sol.q = 2.5;  sol.DP = 0;
x = ones(sol.nel, 1);
sol = init(sol, x);
sol = newt(sol);
plotFigs(sol, x, 0)

% figure;
% eldraw2(sol.ex, sol.ey, [1 2 1]);
% hold on
% eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);
dof = 13;
fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

%% Finite diff
h = 1e-6;

% c = [0.3 0.5 0.2 0.7 0.9]';
% x = repmat(c, sol.nel/5, 1);
x = rand(sol.nel, 1);
sol = init(sol, x);
sol = newt(sol);
[~, dg0, ~, ~] = funcEval(sol, x);

errperc = zeros(sol.nel, 1);
for el = 1:sol.nel
    sol1 = Solver(params);
    sol2 = Solver(params);
    x1 = x;
    x2 = x;
    x1(el) = x1(el) - h;
    x2(el) = x2(el) + h;

    sol1 = init(sol1, x1);
    sol2 = init(sol2, x2);

    sol1 = newt(sol1);
    sol2 = newt(sol2);

    [g1, ~, ~, ~] = funcEval(sol1, x1);
    [g2, ~, ~, ~] = funcEval(sol2, x2);

    dgf = (g2-g1)/2/h;


    fprintf("El: %i \n", el)
    fprintf("  Diff: %.5g \ndg0: %.5g \ndgf: %.5g\n", [dgf-dg0(el), dg0(el), dgf])
    errperc(el) = (dgf-dg0(el))/dg0(el)*1000;
end

patch(sol.ex', sol.ey', errperc);
colormap jet;
colorbar;
%%
N = 10; 
reverse = 0;

epse = zeros(N,1); % eps eff
sige = zeros(N,1); % sig eff
sig = [0;0;0;0];
eps = [0;0;0;0];
deps = [1e-4, 0, 0, 1e-4]'*10;
sigy = sol.sigy0;
Ds = sol.De;
ep = 0;

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
    eps = eps + deps;
    epse(n) = sqrt(eps'*R*eps);
    depsm = deps + [0;0;0;1]*deps(4);
    sigtr = sol.De*depsm + sig;
    if sqrt(sol.sigy0^2*sigtr'*sol.P*sigtr) > sigy
        if sol.DP
            [sig, Dt, Ds, ep] = DPMat(sol, eps, Ds, ep, 1, 1);
        else
            [sig, Dt, sigy] = EPMat(sol, sigtr, sigy, 1, 1);
        end
    else
        sig = sigtr;
        Dt = sol.De;
    end
    sige(n) = sqrt(sol.sigy0^2*sig'*sol.P*sig);
end

figure;
plot([0; epse], [0; sige]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
title("Hill Deformation Model")
grid on;