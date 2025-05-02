clc, clear, close all

% FEM parameters
params.le = 0.1;
params.lx = 2;
params.ly = 1;
params.Vf = 0.3;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e-1;
params.N = 3;
params.disp = [2 -9e-3;
               4 -9e-3;
               6 -9e-3]; % displacement [nodes total-size]

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;
params.v21 = v; params.v32 = v; params.v31 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;

params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy0; %extra terms linHard (sat. stress)
params.xi = 1e-3; %extra terms linHard (sat. exp)

params.rtol = 1e-1;
params.DP = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

% Optimization Parameters
params.re = 2; % Elements in radius
params.filtOn = true;
params.p = 3;
params.q = 2;
params.del = 1e-9;
params.ncon = 1; % Nr of constraints
params.xtol = 1e-4;
params.iterMax = 50;
params.eta = 0.5;
params.beta = 1e-6+0.5;

sol = Solver(params);

%% Optimization
x = ones(sol.nel, 1);
[sol, x] = optimizer(sol, x);

%% Mesh
patch(sol.ex', sol.ey', ones(sol.nel, 1));
colorbar;
colormap jet;

%% Newton-Raphson
x = load("x.mat");
sol = initOpt(sol, x.x);
sol = newt(sol);
plotFigs(sol, x.x, 1, 0)

figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);
dof = 8;
fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

%% Finite diff
h = 1e-6;

c = [0.3 0.5 0.2 0.7 0.9]';
x = repmat(c, sol.nel/5, 1);
% x = 0.8*ones(sol.nel, 1);
sol = initOpt(sol, x);
sol = newt(sol);
[~, dg0, ~, ~] = funcEval(sol, x);

wrong = [];
for el = 1:sol.nel
    sol1 = Solver(params);
    sol2 = Solver(params);
    x1 = x;
    x2 = x;
    x1(el) = x1(el) - h;
    x2(el) = x2(el) + h;

    sol1 = initOpt(sol1, x1);
    sol2 = initOpt(sol2, x2);

    sol1 = newt(sol1);
    sol2 = newt(sol2);

    [g1, ~, ~, ~] = funcEval(sol1, x1);
    [g2, ~, ~, ~] = funcEval(sol2, x2);

    dgf = (g2-g1)/2/h;
    fprintf("El: %i \n", el)
    fprintf("  Diff: %.5g \ndg0: %.5g \ndgf: %.5g\n", [dgf-dg0(el), dg0(el), dgf])
    if (dgf-dg0(el))/dg0(el) > 1e-4
        wrong = [wrong; el];
    end
end

y = ones(sol.nel, 1);
y(wrong) = 0;
patch(sol.ex', sol.ey', y);

