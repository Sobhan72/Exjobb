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

%Material parameters
E = 210e9; v = 0.3; sigma_y0 = 360e6; H = 10e9; K = 160e9;
%beta = linspace(0,50);
sigma_star = [1,1,1,1,0,0];
sigma_star_eff = 

function G = shearMod(E, v, )
G =  E/(2*(1+v));
