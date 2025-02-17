%Hill elastoplastic
clc, clear
sig_y0 = 360e6; 
Fco = 1/(2*sig_y0^2); Gco = 1/(2*sig_y0^2); Hco = 1/(2*sig_y0^2); Lco = 3/(2*sig_y0^2);
P = [Fco+Gco -Fco -Gco 0; -Fco Fco+Hco -Hco 0; -Gco -Hco Gco+Hco 0 ; 0 0 0 2*Lco];

H = 10e9; E = 210e9; v = 0.3; K = E/(3*(1-2*v)); Ge = E/(2*(1+v)); G = Ge; ep = [2 2 1]; % Isotropic elasticity

% E1 = E; E2 = E; E3 = E;
% v12 = v; v13 = v; v23 = v;
% C = [1/E1, -v12/E1, -v13/E1, 0;
%    -v12/E1, 1/E2, -v23/E2, 0;
%    -v13/E1, -v23/E2, 1/E3, 0;
%     0, 0, 0, 1/G12];
% 
% D = inv(C);

D = hooke(2, E, v);

N = 50;
rtol = 1e-4;
eps_eff = zeros(N,1);
sig_eff = zeros(N,1);
eps_inc = [8e-5, 0, 0, 8e-5]';
eps = zeros(4,1);
sig = zeros(4,1);
e_p_eff = 0;
Gp(sig_eff(1), eps_eff(1), e_p_eff, sig_y0, G, Ge, H, rtol, P)

for n = 1:N
    fprintf("Load step: %i \n", n)
    eps = eps + eps_inc;
    eps_eff(n) = strain_eff(eps);
    sig = update_stress(Ge, K, eps_inc) + sig;
    sig_eff(n) = stress_eff(sig);
    
    if sig_eff(n) > sig_y0
        [G, sig_eff(n), e_p_eff, dG] = Gp(sig_eff(n-1), eps_eff(n), e_p_eff, sig_y0, G, Ge, H, rtol, P);
        sig = update_stress(G, K, eps);
        D = Dtan(sig, sig_y0, stress_eff(sig), D, H, P);
    end
end


plot([0; eps_eff], [0; sig_eff]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
grid on;


function [G, sig_eff, e_p_eff, dG] = Gp(sig_eff, eps_eff, e_p_eff, sig_y0, G, Ge, H, rtol, P)
r = sig_eff - 3*G*eps_eff;
iter = 0;
T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 2];
while norm(r) > rtol
    iter = iter + 1;
    Gm = inv(1/Ge*eye(4) + 2*sig_y0^2/sig_eff*e_p_eff*P/T);
    dGdepm = -Gm*(2*sig_y0^2/sig_eff*P/T)*Gm;
    dGdep = dGdepm(1);
    drdep = H - dGdep;
    delta_e_p_eff = -r/drdep;
    e_p_eff = e_p_eff + delta_e_p_eff;
    sig_eff = sig_y0 + H*e_p_eff;
    Gm = inv(1/Ge*eye(4) + 2*sig_y0^2/sig_eff*e_p_eff*P/T);
    G = Gm(1);
    r = sig_eff - 3*G*eps_eff;
    drde = -3*G;
    dG = -drde/drdep;
    fprintf("  iter: %i, r: %4.2g \n", [iter, norm(r)])
end
end

function sig = update_stress(G, K, eps)
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
sig_kk = 3*K*sum(eps(1:3));
s = 2*G*e;
sig = [s(1:3)+sig_kk/3; s(4)];
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