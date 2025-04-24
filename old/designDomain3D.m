function [coord, dof, edof, Ex, Ey, Ez, bc] = designDomain3D(Lx, Ly, Lz, le)
% Generates a recangular design domain with length Lx in the x-direction
% and Ly in the y-direction. Four node elements are used with the side
% length le (same in both x- and y-direction). You can check the mesh by
% running patch(Ex', Ey', 1) in the command window. 
% Note that symmetry is considered along x=0. The symmetry boundary conditions
% are imposed via the matrix bc which contains the degrees of freedom where
% displacements are prescribed in the first column, and the value in the
% seconed column. Extend this matrix with other desirable boundary conditions.

elem_x = Lx/le; elem_y = Ly/le; elem_z = Lz/le; nelm = elem_x*elem_y*elem_z;
nodes_x = elem_x + 1; nodes_y = elem_y + 1; nodes_z = elem_z +1; nnod = nodes_x*nodes_y*nodes_z;
%coord
coord = zeros(nnod, 3);
node  = 1;
for z = 0:nodes_z-1
    for y = 0:nodes_y-1
        for x = 0:nodes_x-1
        coord(node, :) = [x*le y*le z*le];
        node = node + 1;
        end
    end
end
%coord done

%dof
dof = zeros(nnod, 3);
for i = 1:nnod
   dof(i, :) = [i*3-2 i*3-1 i*3]; 
end
%dof done

%enod
enod = zeros(nelm, 9);
enod(:,1) = 1:nelm;
enod(1,2:9) = [1 2 nodes_x+2 nodes_x+1 (nodes_x*nodes_y + 1) (nodes_x*nodes_y + 2) (nodes_x*nodes_y + nodes_x + 2) (nodes_x*nodes_y + nodes_x + 1)];
for i = 2:nelm
   if  mod(i, elem_x*elem_y) == 1
     enod(i, 2:9) = enod(i-1, 2:9) + nodes_x + 2;
   elseif (mod(i, elem_x) == 1)
      enod(i, 2:9) = enod(i-1, 2:9) + 2; 
   else
      enod(i, 2:9) = enod(i-1, 2:9) + 1;
   end
end
%enod done
%Ex, Ey
Ex = zeros(nelm, 8);
Ey = zeros(nelm, 8);
Ez = zeros(nelm, 8);
for i = 1:nelm
   Ex(i, :) = coord(enod(i, 2:9), 1);
   Ey(i, :) = coord(enod(i, 2:9), 2);
   Ez(i, :) = coord(enod(i, 2:9), 3);
end
%Ex, Ey done

%edof
edof = zeros(nelm, 25);
edof(:,1) = 1:nelm;
for i = 1:nelm
   edof(i, 2:25) = [enod(i, 2)*3-2 enod(i, 2)*3-1 enod(i, 2)*3, ...
                  enod(i, 3)*3-2 enod(i, 3)*3-1 enod(i, 3)*3, ...
                  enod(i, 4)*3-2 enod(i, 4)*3-1 enod(i, 4)*3, ...
                  enod(i, 5)*3-2 enod(i, 5)*3-1 enod(i, 5)*3, ...
                  enod(i, 6)*3-2 enod(i, 6)*3-1 enod(i, 6)*3, ...
                  enod(i, 7)*3-2 enod(i, 7)*3-1 enod(i, 7)*3, ...
                  enod(i, 8)*3-2 enod(i, 8)*3-1 enod(i, 8)*3, ...
                  enod(i, 9)*3-2 enod(i, 9)*3-1 enod(i, 9)*3];
end
%edof done

%symmetry bc
bc = [];
for i = 1:nnod
   if (coord(i, 1) == 0.0)
        bc = [bc ; i*3 - 2 0];
   end
end

end

