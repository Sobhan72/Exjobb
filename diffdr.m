E = 210e9; v = 0.3; sig_y0 = 360e6; H = 10e9; K = E/(3*(1-2*v)); Ge = E/(2*(1+v)); G = Ge;

e_p = [2e-3, 0, 0, 2e-3]';
h = [1e-7, 0, 0, 1e-7]';
eps_p_eff = strain_eff(e_p);
sig_eff = sig_y0 + H*eps_p_eff;

s = (2/3)*e_p*sig_eff/eps_p_eff;
e = (1/(2*Ge) + 3*eps_p_eff/(2*sig_eff))*s;
eps_eff = strain_eff(e);


e_p1 = e_p + h;
eps_p_eff1 = strain_eff(e_p1);
sig_eff1 = sig_y0 + H*eps_p_eff1;

s1 = (2/3)*e_p1*sig_eff1/eps_p_eff1;
e1 = (1/(2*Ge) + 3*eps_p_eff1/(2*sig_eff1))*s1;
eps_eff1 = strain_eff(e1);


r = sig_y0 + H*eps_p_eff - 3*Ge/(1+Ge*3*eps_p_eff/(sig_y0 + H*eps_p_eff))*eps_eff;

r1 = sig_y0 + H*eps_p_eff1 - 3*Ge/(1+Ge*3*eps_p_eff1/(sig_y0 + H*eps_p_eff1))*eps_eff1;

dr = H - 3*G - 9*Ge^2*(eps_p_eff*H - sig_eff)/((3*Ge*eps_p_eff + sig_eff)^2)*eps_eff;

(r1-r)/norm(h) -dr



function strain_eff_out = strain_eff(eps)    % Calculate effective stress
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
strain_eff_out = sqrt(2/3*(e'*e + e(4)*e(4)));
end