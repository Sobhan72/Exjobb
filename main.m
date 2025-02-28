clc, clear, close all

% Input parameters
params.le = 1;
params.lx = 1;
params.ly = 1;

params.E = 210e9;
params.v = 0.3;
params.ptype = 2;
params.t = 1;
params.ir = 2;
params.sig_y0 = 360e6;
params.H = 10e9;
params.r2tol = 1e-5;

params.E1 = params.E; params.E2 = params.E; params.E3 = params.E;
params.v12 = params.v; params.v13 = params.v; params.v23 = params.v;
params.v21 = params.v; params.v32 = params.v; params.v31 = params.v;

params.Fco = 1/(2*params.sig_y0^2); 
params.Gco = 1/(2*params.sig_y0^2); 
params.Hco = 1/(2*params.sig_y0^2); 
params.Lco = 3/(2*params.sig_y0^2);

params.N = 10;
params.r1tol = 1e-5;
params.disp = [2 -5e-3]; % displacement [nodes total-size]
sol = Solver(params);

%% Mesh
patch(sol.ex', sol.ey', 1)

%% Newton-Raphson
sol = newt(sol);

figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);

%% FEM test
sol.r1 = ones(sol.ndof, 1);
sol = FEM(sol, params.disp);
figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);

%% Hill model test
N = 50;

epse = zeros(N,1); % eps eff
sige = zeros(N,1); % sig eff
siggp = [0;0;0;0];
epsgp = [0;0;0;0];
deps = [8e-5, 0, 0, 8e-5]';
sigegp = 0;
Dsgp = sol.De;
epgp = 0;

for n = 1:N
    fprintf("Load step: %i \n", n)
    epsgp = epsgp + deps;
    epse(n) = strain_eff(epsgp);
    depsm = deps + [0;0;0;1]*deps(4);
    epsgpm = epsgp + [0;0;0;1]*epsgp(4);
    [siggp, Dtgp, sigegp, Dsgp, epgp] = hill(sol, depsm, epsgpm, siggp, sigegp, Dsgp, epgp);
    sige(n) = sigegp;
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