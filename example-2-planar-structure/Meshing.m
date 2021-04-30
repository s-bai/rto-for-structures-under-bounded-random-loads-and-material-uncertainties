function mesh_info = Meshing(nelx, nely)
    %% Meshing information for sparse matrix (NEW)

    % Input:
    %   nelx, nely: Num. of elements in direction x and y;

    % Output:
    %   mesh_info.nelx: Num. of elements in x direction;
    %   mesh_info.nely: Num. of elements in y direction;
    %   mesh_info.xy: Nodes coordinates
    %   mesh_info.ik & jk: The row and column index vectors
    %   mesh_info.element_node_num: Element node numbers;
    %
    %   Element nodes are arranged as:
    %
    %   i---j
    %   | N |
    %   l---k
    %
    %   Structure is meshed in the order as:
    %
    %       1-----5------9......
    %       | (1) | (4)  |
    %       2-----6-----10......
    %       | (2) | (5)  |
    %       3-----7-----11......
    %       | (3) | (6)  |
    %       4-----8-----12......
    %
    %  x2
    %  |
    %  ---x1  (0,0) Origin is set at node No. 1.

    %% Calculating nodes coordinates
    % The x and y coordinates are stored in the two columns of "xy"
    xy = zeros((nelx + 1) * (nely + 1), 2);
    x_num = 0:nelx;
    y_num = fliplr(0:nely);
    xy(:, 1) = kron(x_num, ones(1, nely + 1))';
    xy(:, 2) = repmat(y_num', nelx + 1, 1);

    %% Calculating node numbers
    node_num = reshape(1:(1 + nelx) * (1 + nely), 1 + nely, 1 + nelx);

    %% Calculating the first node number of each element
    first_node_num_output = reshape(node_num(2:end, 1:end - 1), nelx * nely, 1);
    first_node_num = reshape(2 * node_num(1:end - 1, 1:end - 1) + 1, nelx * nely, 1);

    %% Calculating node numbers of each element
    element_node_num = repmat(first_node_num_output, 1, 4) + repmat([0 nely + 1 nely -1], nelx * nely, 1);

    %% Calculating the first DOF numbers of all the elements (in vector form)
    % first_DOF_vector = reshape(2*node_num(1:end - 1, 1:end - 1) - 1, nelx*nely, 1);

    %% Calculating DOF numbers of all the elements (stored in each row)
    % DOF_matrix = repmat(first_DOF_vector, 1, 8) + repmat([0 1 2*(nely + 1) + [0 1 2 3] 2 3], nelx*nely, 1);
    DOF_matrix = repmat(first_node_num, 1, 8) + repmat([0 1 2 * nely + [2 3 0 1] -2 -1], nelx * nely, 1);

    %% Calculating the row and column index vectors
    % ik = reshape(kron(DOF_matrix, ones(8, 1))', 64*nelx*nely, 1);
    % jk = reshape(kron(DOF_matrix, ones(1, 8))', 64*nelx*nely, 1);

    ik = reshape(kron(DOF_matrix, ones(8, 1))', 64 * nelx * nely, 1);
    jk = reshape(kron(DOF_matrix, ones(1, 8))', 64 * nelx * nely, 1);

    %% Output
    mesh_info.nelx = nelx;
    mesh_info.nely = nely;
    mesh_info.xy = xy;

    mesh_info.ik = ik;
    mesh_info.jk = jk;

    mesh_info.element_node_num = element_node_num;

end % End of Meshing()
