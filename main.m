%von Mises deformation plasticity
clear, clc;
%Making the mesh
L = 0.01;
le = 0.5*L;
lx = 10*L;
ly = 5*L;
[coord, dof, enod, edof, ex, ey, bc] = designDomain(lx, ly, le);
nelm = 

%Material parameters
E = 210e9; v = 0.3; sigma_y0 = 360e6; H = 10e9;
beta = linspace(0,)