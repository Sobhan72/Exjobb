clc, clear, close all

% Input parameters
params.le = 0.05;
params.lx = 1;
params.ly = 1;

params.E = 210e9;
params.v = 0.3;
params.t = 2;
params.ngp = 4;
params.sigy0 = 360e6;
params.H = 10e9;
params.r2tol = 1e-5;
params.DP = 0; % 1 for Deformation plasticity, 0 for Elastoplasticity

params.E1 = params.E; params.E2 = params.E; params.E3 = params.E;
params.v12 = params.v; params.v13 = params.v; params.v23 = params.v;
params.v21 = params.v; params.v32 = params.v; params.v31 = params.v;

params.Fco = 1/(2*params.sigy0^2);
params.Gco = 1/(2*params.sigy0^2);
params.Hco = 1/(2*params.sigy0^2);
params.Lco = 3/(2*params.sigy0^2);

params.N = 100;
params.r1tol = 1e-5;
params.disp = [2 -2e-3;
               4 -2e-3;
               6 -2e-3]; % displacement [nodes total-size]
sol = Solver(params);

%% Mesh
patch(sol.ex', sol.ey', 1)

%% Newton-Raphson
sol = newt(sol);

figure;
eldraw2(sol.ex, sol.ey, [1 2 1]);
hold on
eldisp2(sol.ex, sol.ey, sol.ed, [1 4 1], 10);
dof = sol.ndof-3;
fprintf("Disp DOF %i: %.4g \n", [dof, sol.a(dof)])

%% Model validation
params.DP = 0;
sol = Solver(params);
sol = newt(sol);
dof = sol.ndof/2;
disp1 = sol.a(dof);

params.DP = 1;
sol = Solver(params);
sol = newt(sol);
disp2 = sol.a(dof);
fprintf("\nDisp1: %.5g \nDisp2: %.5g \nDiff: %.5g%% \n", [disp1, disp2, (1-disp2/disp1)*100])

%% Hill material model test
N = 70;
reverse = 1;

epse = zeros(N,1); % eps eff
sige = zeros(N,1); % sig eff
siggp = [0;0;0;0];
epsgp = [0;0;0;0];
deps = [1e-4, 0, 0, 1e-4]';
sigegp = 0;
sigygp = 360e6;
Dsgp = sol.De;
epgp = 0;

for n = 1:N
    if n == 30 && reverse
        deps = -deps;
    elseif n == 40 && reverse
        deps = -deps;
    end
    fprintf("Load step: %i \n", n)
    epsgp = epsgp + deps;
    epse(n) = strain_eff(epsgp);
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