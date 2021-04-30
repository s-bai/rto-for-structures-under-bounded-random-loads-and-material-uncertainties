function [tensor_product_matrix, n_row] = self_tensor_prod(vector, n)
    %% This function calculates vector tensor product of itself
    % Based on "setprod" by Mukhtar Ullah
    % Written by Song Bai

    % Input:
    %   vector: Vector that will be used to calculate the tensor product
    %   n: Order of tensor product

    % Output:
    %   tensor_product_matrix: Tensor product in matrix form
    %   n_row: Rows num. of tensor product, also the number of tensor product nodes.

    temp_grid_cell = cell(1, n);

    for ii = 1:n
        temp_grid_cell{1, ii} = vector;
    end

    [temp_grid_matrix{1:n}] = ndgrid(temp_grid_cell{:});

    for i = n:-1:1
        tp_grid_index_unsorted(:, i) = temp_grid_matrix{i}(:);
    end

    tensor_product_matrix = unique(tp_grid_index_unsorted, 'rows');

    if nargout > 1
        [n_row, ~] = size(tensor_product_matrix);
    end

end % End of self_tensor_prod()
