clc, clear, close all

for n = 1:sol.N
    fprintf("Load step: %i \n", n);
    bcD = sol.disp;
    Nr = 0;
    while norm(sol.r1) > sol.r1tol || Nr == 0
        Nr = Nr + 1;
        sol = FEM(sol, bcD);
        bcD(:, 2) = bcD(:, 2)*0;
        fprintf("  Nr: %i, r1: %4.2g \n", [Nr, norm(sol.r1)]);
    end
    sol.eps = sol.epsi; sol.sig = sol.sigi; sol.ep = sol.epi; sol.Ds = sol.Dsi; sol.sigy = sol.sigyi;
end

gp = 12;
p = 6;
q = 5;
d = 1e-9;
x = 0.5;
h = 1e-8;

k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';

dgam = p*(1-d)*x^(p-1); 
dphi = q*(1-d)*x^(q-1);
th = (dgam*phi-dphi*gam)/(phi)^2;
dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;
dPdx = -2/phi^3*sol.P;
depstdx = sol.eps(gp, :)*dDsdx*sol.P/phi^2*Ds*sol.eps(gp, :)' - sol.eps(gp, :)*Ds*sol.P*2*dphi/phi^3*Ds*sol.eps(gp, :)' + sol.eps(gp, :)*Ds*sol.P/phi^2*dDsdx*sol.eps(gp, :)';
dR2dx = dphi*(sol.sigy0+sol.H*sol.ep(gp))-dphi*sol.sigy0*sqrt(epst)-phi*sol.sigy0/2/sqrt(epst)*depstdx;

x = 0.5 - h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';
r1 = phi*(sol.sigy0 + sol.H*sol.ep(gp)) - phi*sol.sigy0*sqrt(epst);

x = 0.5 + h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;
epst = 1/phi^2*sol.eps(gp, :)*Ds*sol.P*Ds*sol.eps(gp, :)';
r2 = phi*(sol.sigy0 + sol.H*sol.ep(gp)) - phi*sol.sigy0*sqrt(epst);

dr = (r2-r1)/(2*h);

dr - dR2dx




k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
dgam = p*(1-d)*x^(p-1); 
dphi = q*(1-d)*x^(q-1);
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
h = 1e-8;

th = (dgam*phi-dphi*gam)/(phi)^2;

dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;

x=0.5-h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds = gam*sol.De*V*sol.De;

x=0.5+h;
gam = d + (1-d)*x^p;
phi = d + (1-d)*x^q;
V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
Ds2 = gam*sol.De*V*sol.De;

dDs = (Ds2-Ds)/(2*h);
dDsdx-dDs;


tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        dgam = p*(1-d)*x^(p-1); 
        dphi = q*(1-d)*x^(q-1);
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        
        th = (dgam*phi-dphi*gam)/(phi)^2;
        
        dDsdx = dgam*sol.De*V*sol.De - gam*sol.De*V*(th*k0*sol.De*sol.P*sol.De)*V*sol.De;

    ke = ke + B'*dDsdx([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
Kt = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g0t = params.disp(:, 2)'*Kt(params.disp(:, 1),:)*disp2;

x = 0.5 - h;
tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        Ds = gam*sol.De*V*sol.De;

    ke = ke + B'*Ds([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
K1 = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g01 = params.disp(:, 2)'*K1(params.disp(:, 1),:)*disp2;

x = 0.5 + h;
tripK = zeros(64,3);
el = 1;
ke = zeros(8);
for gp = 1:sol.ngp
    [B, J] = NablaB(sol, gp, el);
    ixM = 4*sol.ngp*(el-1) + (gp-1)*4 + 1:4*sol.ngp*(el-1) + gp*4;

        k0 = sol.ep(gp)*sol.sigy0^2/(sol.sigy0+sol.H*sol.ep(gp));
        gam = d + (1-d)*x^p;
        phi = d + (1-d)*x^q;
        V = inv(sol.De + gam/phi*k0*sol.De*sol.P*sol.De);
        Ds = gam*sol.De*V*sol.De;

    ke = ke + B'*Ds([1 2 4],[1 2 4])*B*J*sol.t;
end
[rows, cols] = ndgrid(sol.edof(el, :));
tripK((el-1)*64+1:el*64,:) = [rows(:), cols(:), ke(:)];
K2 = sparse(tripK(:,1), tripK(:,2), tripK(:,3), sol.ndof, sol.ndof);
g02 = params.disp(:, 2)'*K2(params.disp(:, 1),:)*disp2;

dg = (g02-g01)/(2*h);

g0t-dg;

%% Func.
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