function [coord, dof, enod, edof, Ex, Ey, bc] = designDomain(Lx, Ly, le, Wx, Wy)
% Generates a recangular design domain with length Lx in the x-direction
% and Ly in the y-direction. Four node elements are used with the side
% length le (same in both x- and y-direction). You can check the mesh by
% running patch(Ex', Ey', 1) in the command window. 
% Note that symmetry is considered along x=0. The symmetry boundary conditions
% are imposed via the matrix bc which contains the degrees of freedom where
% displacements are prescribed in the first column, and the value in the
% seconed column. Extend this matrix with other desirable boundary conditions.
if isempty(Wx) && isempty(Wy)
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
    edof = zeros(nelm, 8);
    for i = 1:nelm
        edof(i, :) = [enod(i, 2)*2-1 enod(i, 2)*2 enod(i, 3)*2-1 enod(i, 3)*2 enod(i, 4)*2-1 enod(i, 4)*2 enod(i, 5)*2-1 enod(i, 5)*2];
    end
    %edof done

    %symmetry bc
    bc = [];
    for i = 1:nnod
        if (coord(i, 1) == 0.0)
            bc = [bc ; i*2 - 1 0];
        end
    end

else
    % Generates an L-shaped design domain mesh with 4-node elements.
    % Inputs:
    %   Lx - Length of horizontal leg
    %   Ly - Length of vertical leg
    %   le - Element length (square)
    %   Wx - Width of horizontal leg
    %   Wy - Width of vertical leg
    
    % Mesh size
    nx = round(Lx / le);
    ny = round(Ly / le);
    coord = [];
    node_map = zeros(ny+1, nx+1);
    node = 1;

    % Generate coordinates and node map
    for j = 0:ny
        for i = 0:nx
            x = i * le;
            y = j * le;

            % Only include if inside L-beam shape
            inVertical = (x <= Wy) && (y <= Ly);
            inHorizontal = (x >= Wy) && (y <= Wx);

            if inVertical || inHorizontal
                coord(node,:) = [x, y];
                node_map(j+1, i+1) = node;
                node = node + 1;
            else
                node_map(j+1, i+1) = 0;
            end
        end
    end

    nnod = size(coord, 1);

    % DOFs
    dof = zeros(nnod, 2);
    for i = 1:nnod
        dof(i,:) = [2*i - 1, 2*i];
    end

    % Generate elements
    enod = [];
    elem = 1;
    for j = 1:ny
        for i = 1:nx
            n1 = node_map(j, i);
            n2 = node_map(j, i+1);
            n3 = node_map(j+1, i+1);
            n4 = node_map(j+1, i);

            if all([n1 n2 n3 n4] > 0)
                enod(elem,:) = [elem, n1, n2, n3, n4];
                elem = elem + 1;
            end
        end
    end

    nelm = size(enod,1);

    % edof
    edof = zeros(nelm, 8);
    for i = 1:nelm
        n = enod(i,2:5);
        edof(i,:) = reshape([2*n - 1; 2*n], 1, []);
    end

    % Ex, Ey
    Ex = zeros(nelm, 4);
    Ey = zeros(nelm, 4);
    for i = 1:nelm
        Ex(i,:) = coord(enod(i,2:5),1);
        Ey(i,:) = coord(enod(i,2:5),2);
    end

    % Boundary conditions: fix top edge (y = L) and disp on x = L
    L = max(coord(:,2));
    bc = [];
    for i = 1:nnod
        if abs(coord(i,2) - L) < 1e-6
            bc = [bc; 2*i - 1, 0; 2*i, 0];
        end

        if abs(coord(i,1) - L) < 1e-6 && abs(coord(i,2) - Wx/2) <= 2*le
            bc = [bc; 2*i, 1];
        end
        
    end
end
end