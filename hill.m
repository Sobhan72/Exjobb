%Hill elastoplastic
clc, clear, close all
sig_y0 = 360e6; 
F = 1/(2*sig_y0); G = 1/(2*sig_y0); H = 1/(2*sig_y0); L = 3/(2*sig_y0);
P = [F+G -F -G 0; -F F+H -H 0; -G -H G+H 0 ; 0 0 0 2*L];

H = 10e9; E =210e9; v = 0.3; K = E/(3*(1-2*v));
E1 = E; E2 = E; E3 = E;
v12 = v; v13 = v; v23 = v;
G12 = E/(2*(1+v));  %ändra

C = [1/E1, -v12/E1, -v13/E1, 0;
   -v12/E1, 1/E2, -v23/E2, 0;
   -v13/E1, -v23/E2, 1/E3, 0;
    0, 0, 0, 1/G12];

D = inv(C);

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


end

function D_ep_out = D_ep(sig, sig_y0, sig_eff, D, H)
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
A = H + (sig_y0^2/sig_eff)^2*s'*P*D*P*s; %sig_y = sig_eff
D_p = 1/A*(sig_y0^2/sig_eff)^2*D*P*s*s'*P*D; 
D_ep_out = D_el-D_p;
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