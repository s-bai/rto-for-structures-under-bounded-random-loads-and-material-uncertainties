function [TP_nodes, TP_weights] = tp_nodes_generation(n_rand_var, n_tp_grid, GH_nodes, GH_weights)
    %% Generating TENSOR PRODUCT quadrature nodes and weights

    % Input:
    %   n_rand_var: Num. of random variables;
    %   n_tp_grid: Num. of quadrature nodes in each dimension;
    %   GH_nodes: Gauss-Hermite quadrature nodes;
    %   GH_weights: Gauss-Hermite quadrature nodes weights;

    % Output:
    %   TP_nodes: TENSOR PRODUCT quadrature nodes;
    %   TP_weights: TENSOR PRODUCT quadrature nodes weights;

    %% 1-dimensional tensor product base nodes indexes
    tp_base_nodes_idx = 1:n_tp_grid;

    %% Generating tensor product grid nodes indexes
    [tp_nodes_idx, n_tp_nodes] = self_tensor_prod(tp_base_nodes_idx, n_rand_var);

    %% Calculating quadrature nodes coordinates and weights
    TP_nodes = zeros(n_tp_nodes, n_rand_var);
    node_weights = zeros(n_tp_nodes, n_rand_var);

    for ii = 1:n_tp_nodes
        TP_nodes(ii, :) = GH_nodes(tp_nodes_idx(ii, :));
        node_weights(ii, :) = GH_weights(tp_nodes_idx(ii, :));
    end

    %% Change of variables for standard normal distribution
    % Note that the weights here are corrected with '1/(sqrt(pi))^n_rand_var'
    TP_nodes = sqrt(2) * TP_nodes;
    TP_weights = 1 / (sqrt(pi))^n_rand_var * prod(node_weights, 2);

end % End of tp_nodes_generation()
