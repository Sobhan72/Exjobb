clc, clear, close all

% FEM parameters
params.le = 0.001;
params.lx = 0.16; %0.1;
params.ly = 0.16; %0.05;
params.wx = []; %0.04;
params.wy = []; %0.04;
params.loadcase = 4; %1;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e-2;
params.N = 3; % Loadsteps
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
params.re = 3; % Elements in radius
params.filtOn = true; % Filter on
params.pad = true; % Padding on displacement BC (L-beam)

params.p = 1.5;
params.q = 1;
params.del = 1e-9; 
params.dels = 1e-3;
params.rampPQ = [3, 0.1+2/30]; % [end value of p, increment size]

params.beta = 0.1;
params.eta = 0.5;
params.rampB = [2, 10, 1.13]; % [0/1/2, end value, factor size]  (0: off, 1: on, 2: on after simp converges)

params.Vf = 0.3;
params.xtol = 1e-3;
params.ftol = 0.2;
params.iterMax = 1000;

params.stressCon = 0;
params.pnm = 8; % p-norm exponent
params.sigc = 1.05; % Stress constraint factor: sigm = sigy0*sigc
params.stressFree = 0; % Width of area in elements left of right boundary where stress is ignored for L-beam
params.zeroGrad = false; % Manually zero stress gradient in stressFree zone.
params.plasticFree = 0; % Width of area in elements left of right boundary where plasticity is ignored for L-beam
params.mma = [0.5, 10, 0.01]; % initial values [move, lower, upper] 
params.mmaEnd = [0, 0.5, 10, 0.01]; % values after iter [iter, move, lower, upper]    

params.prints = [0,0,0]; % [Load step, R1, R2] 
params.plots = 1;
params.saveName = "";
sol = Solver(params);

%% Optimization
x = 0.4*ones(sol.nel, 1);
[sol, x] = opt(sol, x);
saveData(sol, x, params, "data");

%% Draw Design
sol = Solver(params);
sol.drawDesign(sol, val, x, 1);

%% Draw All Designs
clc, clear, close all
JOB = "1723194";
Solver.drawMultipleDesigns(JOB)

%% Mesh with Padding and Displacement for L-beam
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

%% Mesh with Padding and Displacement for Symmetrical Cantilever 
x = ones(sol.nel, 1);
rho = sol.Z*x + sol.Zp;
patch(sol.ex', sol.ey', rho);
hold on
padx = sol.pc(:, 1)' + [-1  1  1 -1]'*params.le/2;
pady = sol.pc(:, 2)' + [-1 -1  1  1]'*params.le/2;
patch(padx, pady , sol.pc(:, 3))
edge = [0 0; params.lx 0; params.lx params.ly; 0 params.ly; 0 0];
plot(edge(:, 1), edge(:, 2), 'k-', 'LineWidth', 2)
el = ismember(sol.edof, sol.disp(:, 1));
plot(sol.ex(el(:, [2 4 6 8]))', sol.ey(el(:, [2 4 6 8]))', 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
bc = ismember(sol.edof, sol.bcS(:, 1));
plot(sol.ex(bc(:, [2 4 6 8]))', sol.ey(bc(:, [2 4 6 8]))', 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'g')
axis equal

%% Newton-Raphson
% params.plasticFree = 0;
sol = Solver(params);
sol = sol.assignVar(val, sol);
sol.p = 4; sol.q = 3.5;
sol = init(sol, x);
sol = newt(sol);
[g0, dg0, gc, dgc] = funcEval(sol, x);
sol.drawDesign(sol, val, x, 1);
sol.plotGrads(dg0, dgc);

fprintf("g0; %.4g \n",  -sol.a(sol.disp(:, 1))'*sol.R1(sol.disp(:, 1)));
% figure;
% eldraw2(sol.ex, sol.ey, [1 2 1]);
% hold on
% eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 5);
% dof = 16;
% fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)]);

%% Finite Difference
h = 5e-6;
% params.plasticFree = 0;

% c = [0.3 0.5 0.01 0.7]';
% x = repmat(c, sol.nel/4, 1);
% x = rand(sol.nel, 1);
sol = Solver(params);
sol = sol.assignVar(val, sol);
sol.p = 4; sol.q = 3.5;
sol0 = sol;
sol0 = init(sol0, x);
sol0 = newt(sol0);
[~, dg, ~, dgc] = funcEval(sol0, x);

wrong = [];
% [~, els] = maxk(dgc(2, :), 10);
els = find(all((sol.ex>0.0365 & sol.ex<0.0395) & (sol.ey<0.0995 & sol.ey>0.0975), 2))';
for el = els
    sol1 = sol;
    sol2 = sol;
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

    dgcf = (gc2(1)-gc1(1))/2/h;
    dgf = (g2-g1)/2/h;

    fprintf("El: %i x: %.5g \n", [el, x(el)])
    fprintf("  Diff: %.5g \ndgc2: %.5g \ndgcf: %.5g\n", [dgcf-dgc(1, el), dgc(1, el), dgcf])
    fprintf("  Diff: %.5g \ndg: %.5g \ndgf: %.5g\n\n", [dgf-dg(el), dg(el), dgf])
    % if (abs((dgcf-dgc(2, el))/dgc(2, el))) > 1e-3 && abs(dgcf) > 5e-6
    %     wrong = [wrong; el];
    % end
end

y = ones(sol.nel, 1);
y(wrong) = 0;
patch(sol.ex', sol.ey', y);