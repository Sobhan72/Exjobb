clc,clear,close all

function D = Dloop(obj, Ds, dDsdep, depdeps, eps)
indx = [1 1 1 1;
    1 1 2 2;
    1 1 3 3;
    1 1 1 2;
    2 2 2 2;
    2 2 3 3;
    2 2 1 2;
    3 3 3 3;
    3 3 1 2;
    1 2 1 2];

Dti = [];

eps = [eps(1) eps(4)/2 0;
    eps(4)/2 eps(2) 0;
    0 0 eps(3)];

depdeps = [depdeps(1) depdeps(4) 0;
    depdeps(4) depdeps(2) 0;
    0 0 depdeps(3)];

dDsdep = tensor(obj, dDsdep);

for ii = indx'
    Dtii = 0;
    for k = 1:3
        for l = 1:3
            Dtii = Dtii + dDsdep(ii(1), ii(2), k, l)*depdeps(ii(3), ii(4))*eps(k,l);
        end
    end
    Dti = [Dti, Dtii];
end

D = Ds + [Dti(1), Dti(2), Dti(3), Dti(4);
    Dti(2), Dti(5), Dti(6), Dti(7);
    Dti(3), Dti(6), Dti(8), Dti(9);
    Dti(4), Dti(7), Dti(9), Dti(10)];
end

function t = tensor(obj, matrix)
t = zeros(3,3,3,3);
t(1,1,1,1) = matrix(1,1);
t(1,1,1,2) = matrix(1,4);
t(1,1,2,1) = matrix(1,4);
t(1,1,2,2) = matrix(1,2);
t(1,1,3,3) = matrix(1,3);

t(1,2,1,1) = matrix(4,1);
t(1,2,1,2) = matrix(4,4);
t(1,2,2,1) = matrix(4,4);
t(1,2,2,2) = matrix(4,2);
t(1,2,3,3) = matrix(4,3);

t(2,1,1,1) = matrix(4,1);
t(2,1,1,2) = matrix(4,4);
t(2,1,2,1) = matrix(4,4);
t(2,1,2,2) = matrix(4,2);
t(2,1,3,3) = matrix(4,3);

t(2,2,1,1) = matrix(2,1);
t(2,2,1,2) = matrix(2,4);
t(2,2,2,1) = matrix(2,4);
t(2,2,2,2) = matrix(2,2);
t(2,2,3,3) = matrix(2,3);

t(3,3,1,1) = matrix(3,1);
t(3,3,1,2) = matrix(3,4);
t(3,3,2,1) = matrix(3,4);
t(3,3,2,2) = matrix(3,2);
t(3,3,3,3) = matrix(3,3);
end

function Dt = newDt(obj, sige, Ds, eps, delta_ep) %Not correct
sig = Ds*eps;
dfdsig = obj.sig_y0^2/sige*obj.P*sig;
df2dsig2 = obj.sig_y0^2/sige*obj.P*(eye(4)-obj.sig_y0^2/sige^2*sig*sig'*obj.P);
Da = inv(obj.C + delta_ep*df2dsig2);
A = dfdsig'*Da*dfdsig + obj.H;
Dt = Da - 1/A*Da*(dfdsig*dfdsig')*Da;
end

function sig = update_stress(obj, eps)
e = [eps(1:3)-mean(eps(1:3)); 2*eps(4)];
sigkk = 3*obj.Kmod*sum(eps(1:3));
s = obj.G*e;
sig = [s(1:3)+sigkk/3; s(4)];
end

function sig = update_stress_el(obj, eps)
e = [eps(1:3)-mean(eps(1:3)); eps(4)];
sigkk = 3*obj.Kmod*sum(eps(1:3));
s = 2*obj.Ge*e;
sig = [s(1:3)+sigkk/3; s(4)];
end

function stress_eff_out = stress_eff(~, sig)    % Calculate effective strain
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
J2 = 0.5*(s'*s + s(4)*s(4));
stress_eff_out = sqrt(J2*3);
end

function D_ep = Dtan(obj, sig, sig_eff)   %%FEL%%
s = [sig(1:3)-mean(sig(1:3)); sig(4)];
A = obj.H + (obj.sig_y0^2/sig_eff)^2*s'*obj.P*D*obj.P*s; %sig_y = sig_eff
D_p = 1/A*(obj.sig_y0^2/sig_eff)^2*D*obj.P*s*s'*obj.P*D;
D_ep = D-D_p;
end

function [bc, disp] = tempBc(~, ly, le, ndof)
nR = ly/le + 1;
bc = [[ndof/nR*(0:nR-1)'+1; ndof/nR*(0:nR-1)'+2] zeros(2*nR,1)];
disp = [ndof/nR*(1:nR)'-1, ones(nR,1)*3e-3];
end

function Dtf = FDM(obj, eps, sige, Ds, ep)
delta = 1e-8;

Dtf = zeros(4);
for i = 1:4
    deps = [0; 0; 0; 0];
    deps(i) = delta;
    eps2 = eps - deps;
    eps3 = eps + deps;
    [~, ~, Ds2, ~] = DMat(obj, eps2, sige, Ds, ep);
    [~, ~, Ds3, ~] = DMat(obj, eps3, sige, Ds, ep);
    sig2 = Ds2*eps2;
    sig3 = Ds3*eps3;
    Dtf(:, i) = (sig3-sig2)/2/delta;
end
end