clc, clear, close all

% Input parameters
params.le = 1;
params.lx = 4;
params.ly = 4;

params.E = 210e9;
params.v = 0.3;
params.eparm = [2 1 2];
params.sig_y0 = 360e6;
params.H = 10e9;

params.rtol = 1e-4;

params.Fco = 1/(2*params.sig_y0^2); 
params.Gco = 1/(2*params.sig_y0^2); 
params.Hco = 1/(2*params.sig_y0^2); 
params.Lco = 3/(2*params.sig_y0^2);


sol = Solver(params);
% patch(sol.ex', sol.ey', 1)



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
    [sol, siggp] = hill(sol, deps, epsgp, siggp);
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