clc, clear, close all

% FEM parameters
params.le = 0.001;
params.lx = 0.1; %0.1;
params.ly = 0.1; %0.05;
params.wx = 0.04; %[];
params.wy = 0.04; %[];
params.loadcase = 3; %1;

params.Vf = 0.35;
params.t = 1;
params.ngp = 4;

params.stressCon = 1;
params.pnm = 8; % p-norm exponent
params.sigc = 1; % Stress constraint factor: sigm = sigy0*sigc

params.R1tol = 1e-2; 
params.N = 4; % Loadsteps
params.disp = -1e-3; % Total displacement 

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;
 
params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; % Extra terms linHard (sat. stress)
params.xi = 1e3; % Extra terms linHard (sat. exp)

params.rtol = 1e-4;
params.PT = 1; % 0 for Incremental plasticity, 1 for Deformation plasticity

% Optimization Parameters
params.re = 3; % Elements in radius
params.filtOn = true;
params.p = 1.5;
params.q = 1; 
params.eta = 0.5;
params.beta = 1;
params.rampB = 1; % 0: off, 1: B*1.1, 2: B + 1
params.rampPQ = true;
params.del = 1e-9; 
params.dels = 1e-3;
params.xtol = 1e-5;
params.iterMax = 500;

params.print = [0,0,0]; % [Load step, R1, R2] 
params.saveName = "";
sol = Solver(params);

%% Optimization
x = ones(sol.nel, 1);
[sol, x] = opt(sol, x);
saveData(sol, x, params, "data");

%% Draw Design
% clc, clear, close all
% load("data\OptDesign_case=1_x=0.8_disp=1e-3.mat")
% load("DesignFailcase=1.mat");

sol = Solver(params);
sol = sol.assignVar(val, sol);
sol.beta = 10; sol.p = 3; sol.q = 2.5;
x = sol.he(sol.Z*x);
sol.phi = sol.dels + (1-sol.dels)*x.^sol.q;
plotFigs(sol, x, 0);
% plotFigs(sol, x, 1);

%% Mesh
patch(sol.ex', sol.ey', ones(sol.nel, 1));
axis equal

%% Newton-Raphson
x = ones(sol.nel, 1);
% load("DesignFailcase=1.mat");
% sol = sol.assignVar(val, sol);
% sol = Solver(params);
% sol.beta = 10; sol.p = 3; sol.q = 3; sol.t = 1; 
% sol.beta = val.beta; sol.p = val.p; sol.q = val.q;
sol = init(sol, x);
sol = newt(sol);
plotFigs(sol, x, 0)

fprintf("g0; %.4g \n",  -sol.a(sol.pdof)'*sol.R1(sol.pdof));
figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 1);
% dof = 16;
% fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

%% Finite diff
h = 1e-6;

% c = [0.3 0.5 0.2 0.7]';
% x = repmat(c, sol.nel/4, 1);
% x = 0.5*ones(sol.nel, 1);
load("x.mat");
sol = init(sol, x);
sol = newt(sol);
[~, ~, ~, dgc] = funcEval(sol, x);

wrong = [];
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

    [~, ~, gc1, ~] = funcEval(sol1, x1);
    [~, ~, gc2, ~] = funcEval(sol2, x2);

    dgf = (gc2(2)-gc1(2))/2/h;


    fprintf("El: %i \n", el)
    fprintf("  Diff: %.5g \ndg0: %.5g \ndgf: %.5g\n", [dgf-dgc(2, el), dgc(2, el), dgf])
    if (abs((dgf-dgc(2, el))/dgc(2, el))) > 1e-4
        wrong = [wrong; el];
    end
end

y = ones(sol.nel, 1);
y(wrong) = 0;
patch(sol.ex', sol.ey', y);