%von Mises deformation plasticity
clear, clc;
%Making the mesh
le = 0.003;
lx = 0.03;
ly = 0.03;
lz = le;
[coord, dof, edof, ex, ey, ez, ~] = designDomain(lx, ly, lz, le);
nel = (lx/le)*(ly/le)*(lz/le);
nnod = size(coord, 1);
ndof = nnod*3;
bc = [];
for i = 1:nnod
   if (coord(i, 1) == 0.0)
        bc = [bc ; i*3 - 2 0 ; i*3 - 1 0 ; i*3 0 ];
   end
end

fload = zeros(ndof, 1);
fload2 = zeros(ndof, 1);
nn = 0;
F = 300;
for i = 1:nnod
    if coord(i,1) == lx && abs(coord(i,2) - ly/2) <= 0.005
        nn = nn + 1;
        fload(3*i-1) = -1;
        fload2(3*i-1) = -2/3;
        fload2(3*i-2) = 2/3;
    end
end
fload = fload*F/nn;
fload2 = fload2*F/nn;

K = zeros(ndof);
for el = 1:nel
    Ke = soli8e(ex,ey,ez,ep,D);
    indx = edof(:, 2:end);
    K(indx, indx) = K(indx, indx) + Ke;
end

a = solveq(K, fload, bc);
ed = extract_ed(edof, a);



%%
clc, clear, close all
%Material parameters
E = 210e9; v = 0.3; sig_y0 = 360e6; H = 10e9; K = E/(3*(1-2*v)); Ge = E/(2*(1+v)); G = Ge; ep = [2 2 1];
D = hooke(2, E, v);

N = 50;
rtol = 1e-4;
eps_eff = zeros(N,1);
sig_eff = zeros(N,1);
eps_inc = [8e-5, 0, 0, 8e-5]';
eps = zeros(4,1);
sig = zeros(4,1);
e_p_eff = 0;

for n = 1:N
    fprintf("Load step: %i \n", n)
    eps = eps + eps_inc;
    eps_eff(n) = strain_eff(eps);  
    sig = update_stress(Ge, K, eps_inc) + sig;
    sig_eff(n) = stress_eff(sig);

    if sig_eff(n) > sig_y0
        [G, sig_eff(n), e_p_eff, dG] = Gp(sig_eff(n-1), eps_eff(n), e_p_eff, sig_y0, G, Ge, H, rtol);
        sig = update_stress(G, K, eps);
        D = Dtan(K, G, dG, eps_eff(n), eps);
    end
end

plot([0; eps_eff], [0; sig_eff]/1e6, 'LineWidth', 2);
xlabel('$\epsilon_{eff}$', 'Interpreter', 'latex'); 
ylabel('$\sigma_{eff}$ (MPa)', 'Interpreter', 'latex');
grid on;

function D = Dtan(K, G, dG, eps_eff, eps)
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
Lam = 1/3*[2, -1, -1, 0; -1, 2, -1, 0; -1, -1, 2, 0; 0, 0, 0, 3/2];
llam = [1, 1, 1, 0; 1, 1, 1, 0; 1, 1, 1, 0; 0, 0, 0, 0];
D = K*llam + 2*G*Lam + 4/3/eps_eff*dG*(e*e')*Lam;
end

function [G, sig_eff, e_p_eff, dG] = Gp(sig_eff, eps_eff, e_p_eff, sig_y0, G, Ge, H, rtol)
r = sig_eff - 3*G*eps_eff;
iter = 0;
while norm(r) > rtol
    iter = iter + 1;
    dGdep = 9*Ge^2*(e_p_eff*H - sig_eff)/(3*Ge*e_p_eff + sig_eff)^2*eps_eff;
    drdep = H - dGdep;
    delta_e_p_eff = -r/drdep;
    e_p_eff = e_p_eff + delta_e_p_eff;
    sig_eff = sig_y0 + H*e_p_eff;
    G = Ge/(1+Ge*3*e_p_eff/sig_eff);
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

function stress_eff_out = stress_eff(sig)    % Calculate effective stress
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
J2 = 0.5*(s'*s + s(4)*s(4));
stress_eff_out = sqrt(J2*3);
end

function strain_eff_out = strain_eff(eps)    % Calculate effective stress
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
strain_eff_out = sqrt(2/3*(e'*e + e(4)*e(4)));
end