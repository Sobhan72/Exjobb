clc, clear, close all

% FEM parameters
params.le = 0.001;
params.lx = 0.1; %0.1;
params.ly = 0.1; %0.05;
params.wx = 0.04; %[];
params.wy = 0.04; %[];
params.loadcase = 3; %1;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e-2;
params.N = 5; % Loadsteps
params.disp = -1.6e-3; % Total displacement 

% Material
E = 210e9;
v = 0.4;
params.E1 = E; params.E2 = 0.7*E; params.E3 = 0.5*E;
params.v12 = v; params.v13 = 0.5*v; params.v23 = 0.5*v;

params.sigy01 = 360e6;
params.sigy02 = 300e6;
params.sigy03 = 250e6;
 
params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; % Extra terms linHard (sat. stress)
params.xi = 1e3; % Extra terms linHard (sat. exp)

params.rtol = 1e-4;
params.PT = 1; % 0 for Incremental plasticity, 1 for Deformation plasticity

% Optimization Parameters
params.re = 6; % Elements in radius
params.filtOn = true; % Filter on
params.pad = true; % Padding on displacement BC (L-beam)

params.p = 1.5;
params.q = 1; 
params.del = 1e-9; 
params.dels = 1e-3;
params.rampPQ = [4, 0.1+2/30]; % [end value of p, increment size]

params.beta = 0.1;
params.eta = 0.5;
params.rampB = [2, 10, 1.13]; % [0/1/2, end value, factor size]  (0: off, 1: on, 2: on after simp converges)

params.Vf = 0.3;
params.xtol = 1e-3;
params.iterMax = 1250;

params.stressCon = 1;
params.pnm = 8; % p-norm exponent
params.sigc = 1.15; % Stress constraint factor: sigm = sigy0*sigc
params.stressFree = 15; % Width of area in elements left of right boundary where stress is ignored for L-beam
params.mma = [0.1, 10, 0.01]; % initial values [move, lower, upper] 
params.mmaEnd = [350, 0.05, 0.1, 0.001]; % values after iter [iter, move, lower, upper]    

params.print = [0,0,0]; % [Load step, R1, R2] 
params.plots = 1;
params.saveName = "";
sol = Solver(params);

%% Optimization
x = 0.35*ones(sol.nel, 1);
[sol, x] = opt(sol, x);
saveData(sol, x, params, "data");

%% Draw Design
sol = Solver(params);
sol.drawDesign(sol, val, x, 0);

%% Draw All Designs
clc, clear, close all
JOB = "1657973";
Solver.drawMultipleDesigns(JOB)

%% Mesh with Padding and Displacement
x = ones(sol.nel, 1);
rho = sol.Z*x + sol.Zp;
el = ismember(sol.edof, sol.disp(:, 1));
padx = sol.pc(:, 1)' + [-1  1  1 -1]'*params.le/2;
pady = sol.pc(:, 2)' + [-1 -1  1  1]'*params.le/2;
edge = [0 0; params.lx 0; params.lx params.wx; params.wy params.wx;
        params.wy params.ly; 0 params.ly; 0 0];
patch(sol.ex', sol.ey', rho);
hold on
patch(padx, pady , sol.pc(:, 3))
plot(edge(:, 1), edge(:, 2), 'k-', 'LineWidth', 2)
plot(sol.ex(el(:, [2 4 6 8]))', sol.ey(el(:, [2 4 6 8]))', 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
axis equal

%% Newton-Raphson
% load("DesignFailcase=1.mat");
% sol = sol.assignVar(val, sol);
sol = Solver(params);
sol.beta = 0.01; sol.p = 4; sol.q = 3.5; 
% sol.beta = val.beta; sol.p = val.p; sol.q = val.q;
sol = init(sol, x);
sol = newt(sol);
% plotFigs(sol, x, 0)

fprintf("g0; %.4g \n",  -sol.a(sol.pdof)'*sol.R1(sol.pdof));
figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 5);
% dof = 16;
% fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

%% Finite Difference
h = 5e-6;

c = [0.3 0.5 0.01 0.7]';
x = repmat(c, sol.nel/4, 1);
% x = rand(sol.nel, 1);
% load("x.mat");
sol = init(sol, x);
sol = newt(sol);
[~, dg, ~, dgc] = funcEval(sol, x);

wrong = [];
[~, els] = maxk(dgc(2, :), 10);
for el = els
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

    [g1, ~, gc1, ~] = funcEval(sol1, x1);
    [g2, ~, gc2, ~] = funcEval(sol2, x2);

    dgcf = (gc2(2)-gc1(2))/2/h;
    dgf = (g2-g1)/2/h;

    fprintf("El: %i \n", el)
    fprintf("  Diff: %.5g \ndgc2: %.5g \ndgcf: %.5g\n", [dgcf-dgc(2, el), dgc(2, el), dgcf])
    %fprintf("  Diff: %.5g \ndg: %.5g \ndgf: %.5g\n", [dgf-dg(el), dg(el), dgf])
    if (abs((dgcf-dgc(2, el))/dgc(2, el))) > 1e-3 && abs(dgcf) > 5e-6
        wrong = [wrong; el];
    end
end

y = ones(sol.nel, 1);
y(wrong) = 0;
patch(sol.ex', sol.ey', y);