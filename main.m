%von Mises deformation plasticity
clear, clc;
%Making the mesh
L = 0.01;
le = 0.5*L;
lx = 10*L;
ly = 5*L;
[coord, dof, enod, edof, ex, ey, bc] = designDomain(lx, ly, le);
nelm = (lx/le)*(ly/le);
%patch(ex', ey', 1)

%%
clc, clear, close all
%Material parameters
E = 210e9; v = 0.3; sig_y0 = 360e6; H = 10e9; K = E/(3*(1-2*v)); Ge = E/(2*(1+v)); G = Ge;
D_el = hooke(2, E, v);

N = 50;
rtol = 1e-5;
eps_eff = zeros(N,1);
sig_eff = zeros(N,1);
eps_inc = [1e-4, 0, 0, 1e-4]';
eps = zeros(4,1);
sig_1 = zeros(4,1);
e_p_eff = 0;

for n = 1:N
    eps = eps + eps_inc;
    eps_eff(n) = strain_eff(eps);  
    sig_t = D_el*eps_inc + sig_1;
    sig_t_eff = stress_eff(sig_t);

    if sig_t_eff > sig_y0
        r = sig_eff(n) - 3*G*eps_eff(n);
        sig_eff_iter = sig_eff(n-1);
        while norm(r) > rtol
            dr = H - 3*G - 9*Ge^2*(e_p_eff*H - sig_eff_iter)/(3*Ge*e_p_eff + sig_eff_iter)^2*eps_eff(n);
            delta_e_p_eff = -r/dr;
            e_p_eff = e_p_eff + delta_e_p_eff;
            sig_eff_iter = sig_y0 + H*e_p_eff;
            G = Ge/(1+Ge*3*e_p_eff/sig_eff_iter);
            r = sig_eff_iter - 3*G*eps_eff(n);
        end
        sig_eff(n) = sig_eff_iter;
        sig = update_stress(G, K, eps);
    else
        sig = D_el*eps;
        sig_eff(n) = stress_eff(sig);
        eps_eff(n) = strain_eff(eps);
    end
    sig_1 = sig;
end

plot(eps_eff, sig_eff)

function sig = update_stress(G, K, eps)
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
sig_kk = 3*K*sum(eps(1:3));
s = 2*G*e;
sig = [s(1:3)+sig_kk/3; s(4)];
end

function stress_eff_out = stress_eff(sig)    % Calculate effective stress
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
J2 = 0.5*(s'*s + s(4)*s(4));
stress_eff_out = sqrt(J2*3);
end

function strain_eff_out = strain_eff(eps)    % Calculate effective stress
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
strain_eff_out = sqrt(2/3*e'*e + e(4)*e(4));
end




% sigma_star = [1,1,1,1,0,0]';
% sigma_star_eff = stress_eff(sig_star);
% sigma_star_kk = sum(sigma_star(1:3));
% sig = Beta(eps, K, sigma_star_kk)*sigma_star;
% function b = beta(eps, K, sig_star_kk)
% sig_kk = sum(eps(1:3))*3*K;
% b = sig_kk/sig_star_kk;
% end