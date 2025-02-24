clc, clear, close all

% Input parameters
params.le = 1e-1;
params.lx = 4;
params.ly = 4;

params.E = 210e9;
params.v = 0.3;
params.epm = [2 1 2];
params.sig_y0 = 360e6;
params.H = 10e9;

params.E1 = params.E; params.E2 = params.E; params.E3 = params.E;
params.v12 = params.v; params.v13 = params.v; params.v23 = params.v;
params.v21 = params.v; params.v32 = params.v; params.v31 = params.v;

params.rtol = 1e-4;

params.Fco = 1/(2*params.sig_y0^2); 
params.Gco = 1/(2*params.sig_y0^2); 
params.Hco = 1/(2*params.sig_y0^2); 
params.Lco = 3/(2*params.sig_y0^2);


sol = Solver(params);
% patch(sol.ex', sol.ey', 1)

%% FEM test
bcD = [2 -1e-4];
sol.res = ones(sol.ndof, 1);
sol = FEM(sol, sol.De, bcD);
figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1]);

%% Hill model test
N = 50;

epse = zeros(N,1); % eps eff
sige = zeros(N,1); % sig eff
siggp = [0;0;0;0];
epsgp = [0;0;0;0];
deps = [8e-5, 0, 0, 8e-5]';

for n = 1:N
    fprintf("Load step: %i \n", n)
    epsgp = epsgp + deps;
    epse(n) = strain_eff(epsgp);
    depsm = deps + [0;0;0;1]*deps(4);
    epsgpm = epsgp + [0;0;0;1]*epsgp(4);
    [sol, siggp] = hill(sol, depsm, epsgpm, siggp);
    sige(n) = sol.sig_eff;
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