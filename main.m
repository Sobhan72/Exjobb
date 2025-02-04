%von Mises deformation plasticity
clear, clc;
%Making the mesh
L = 0.01;
le = 0.1*L;
[coord, dof, enod, edof, ex, ey, bc] = designDomain(6*L, 6*L, le); 
patch(ex',ey',1)
