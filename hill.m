%Hill elastoplastic
clc, clear
sig_y0 = 360e6; 
Fco = 1/(2*sig_y0^2); Gco = 1/(2*sig_y0^2); Hco = 1/(2*sig_y0^2); Lco = 3/(2*sig_y0^2);
P = [Fco+Gco -Fco -Gco 0; -Fco Fco+Hco -Hco 0; -Gco -Hco Gco+Hco 0 ; 0 0 0 2*Lco];

H = 10e9; E = 210e9; v = 0.3; K = E/(3*(1-2*v)); Ge = E/(2*(1+v)); mp = [2 2 1]; % Isotropic elasticity

E1 = E; E2 = E; E3 = E;
v12 = v; v13 = v; v23 = v; v21 = v; v32 = v; v31 = v;
G12 = Ge; 
C = [1/E1, -v21/E2, -v31/E3, 0;
   -v12/E1, 1/E2, -v32/E3, 0;
   -v13/E1, -v23/E2, 1/E3, 0;
    0, 0, 0, 1/G12];

De = inv(C); %elastic D matrix
Dp = De; 

N = 50;
rtol = 1e-4;
eps_eff = zeros(N,1);
sig_eff = zeros(N,1);
eps_inc = [8e-5, 0, 0, 8e-5]';
eps = zeros(4,1);
sig = zeros(4,1);
ep = 0; % effective plastic strain

%FDM
% eps = eps + eps_inc;
% e = [eps(1:3)-mean(eps(1:3)); 2*eps(4)];
% h = 1e-9;
% ep1 = 0.001;
% sig_eff = sig_y0 + H*ep1;
% G_inv = 1/2/Ge*T + sig_y0^2/sig_eff*ep1*P;
% G = inv(G_inv);
% dGdep = -G*P*G*(sig_y0^2*(sig_eff-ep1*H)/sig_eff^2);
% detdG = 2*P*G*e*e';
% et1 = e'*G*P*G*e;
% detdep = trace(detdG*dGdep);
% 
% ep1 = ep1 + h;
% sig_eff = sig_y0 + H*ep1;
% G_inv = 1/2/Ge*T + sig_y0^2/sig_eff*ep1*P;
% G = inv(G_inv);
% dGdep = -G*P*G*(sig_y0^2*(sig_eff-ep1*H)/sig_eff^2);
% detdG = 2*P*G*e*e';
% et2 = e'*G*P*G*e;
% detdepf = (et2 - et1)/h;


for n = 1:N
    fprintf("Load step: %i \n", n)
    eps = eps + eps_inc;
    epsm = [eps(1:3); 2*eps(4)];
    eps_eff(n) = strain_eff(eps);
    sig = update_stress_el(De, eps_inc) + sig;
    sig_eff(n) = stress_eff(sig);
    
    if sig_eff(n) > sig_y0
        [Dp, sig_eff(n), ep] = Ds(sig_eff(n-1), epsm, ep, sig_y0, C, Dp, H, rtol, P);
        sig = update_stress(Dp, epsm);
        D = Dtan(sig, sig_y0, stress_eff(sig), De, H, P);
    end
end

figure;
plot([0; eps_eff], [0; sig_eff]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
title("Hill Deformation Model")
grid on;


function [Dp, sig_eff, ep] = Ds(sig_eff, eps, ep, sig_y0, C, Dp, H, rtol, P)
epst = eps'*Dp*P*Dp*eps;
r = sig_eff - sig_y0*sqrt(epst);
iter = 0;
while norm(r) > rtol
    iter = iter + 1;
    Dp = inv(C + sig_y0^2/sig_eff*ep*P);
    dDdep = -Dp*P*Dp*(sig_y0^2*(sig_eff-ep*H)/sig_eff^2);
    detdD = 2*P*Dp*eps*eps';
    epst = eps'*Dp*P*Dp*eps;
    drdep = H - sig_y0/(2*sqrt(epst))*trace(detdD*dDdep);
    delta_ep = -r/drdep;
    ep = ep + delta_ep;
    sig_eff = sig_y0 + H*ep;
    Dp = inv(C + sig_y0^2/sig_eff*ep*P);
    epst = eps'*Dp*P*Dp*eps;
    r = sig_eff - sig_y0*sqrt(epst);
    fprintf("  iter: %i, r: %4.2g \n", [iter, norm(r)])
end
end

function sig = update_stress(Dp, epsm)
sig = Dp*epsm;
end

function sig = update_stress_el(De, eps)
eps = [eps(1:3); 2*eps(4)];
sig = De*eps;
end

function D_ep = Dtan(sig, sig_y0, sig_eff, D, H, P)
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
A = H + (sig_y0^2/sig_eff)^2*s'*P*D*P*s; %sig_y = sig_eff
D_p = 1/A*(sig_y0^2/sig_eff)^2*D*P*s*s'*P*D; 
D_ep = D-D_p;
end

function stress_eff_out = stress_eff(sig)    % Calculate effective strain
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
J2 = 0.5*(s'*s + s(4)*s(4));
stress_eff_out = sqrt(J2*3);
end

function strain_eff_out = strain_eff(eps)    % Calculate effective stress
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
strain_eff_out = sqrt(2/3*(e'*e + e(4)*e(4)));
end