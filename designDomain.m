function [coord, dof, enod, edof, Ex, Ey, bc] = designDomain(Lx, Ly, le)
% Generates a recangular design domain with length Lx in the x-direction
% and Ly in the y-direction. Four node elements are used with the side
% length le (same in both x- and y-direction). You can check the mesh by
% running patch(Ex', Ey', 1) in the command window. 
% Note that symmetry is considered along x=0. The symmetry boundary conditions
% are imposed via the matrix bc which contains the degrees of freedom where
% displacements are prescribed in the first column, and the value in the
% seconed column. Extend this matrix with other desirable boundary conditions.

elem_x = Lx/le; elem_y = Ly/le; nelm = elem_x*elem_y;
nodes_x = elem_x + 1; nodes_y = elem_y + 1; nnod = nodes_x*nodes_y;
%coord
coord = zeros(nnod, 2);
node  = 1;
for y = 0:nodes_y-1
    for x = 0:nodes_x-1
       coord(node, :) = [x*le y*le];
       node = node + 1;
    end
end
%coord done

%dof
dof = zeros(nnod, 2);
for i = 1:nnod
   dof(i, :) = [i*2-1 i*2]; 
end
%dof done

%enod
enod = zeros(nelm, 5);
enod(:,1) = 1:nelm;
enod(1,2:5) = [1 2 nodes_x+2 nodes_x+1];
for i = 2:nelm
   if (mod(i, elem_x) == 1)
      enod(i, 2:5) = enod(i-1, 2:5) + 2; 
   else
      enod(i, 2:5) = enod(i-1, 2:5) + 1;
   end
end
%enod done
%Ex, Ey
Ex = zeros(nelm, 4);
Ey = zeros(nelm, 4);
for i = 1:nelm
   Ex(i, :) = coord(enod(i, 2:5), 1);
   Ey(i, :) = coord(enod(i, 2:5), 2);
end
%Ex, Ey done

%edof
edof = zeros(nelm, 9);
edof(:,1) = 1:nelm;
for i = 1:nelm
   edof(i, 2:9) = [enod(i, 2)*2-1 enod(i, 2)*2 enod(i, 3)*2-1 enod(i, 3)*2 enod(i, 4)*2-1 enod(i, 4)*2 enod(i, 5)*2-1 enod(i, 5)*2]; 
end
%edof done

%symmetry bc
bc = [];
for i = 1:nnod
   if (coord(i, 1) == 0.0)
        bc = [bc ; i*2 - 1 0];
   end
end

end