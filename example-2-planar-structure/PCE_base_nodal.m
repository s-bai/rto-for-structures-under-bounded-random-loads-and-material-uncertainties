function PCE_base_nodal_value = PCE_base_nodal(Q_nodes, pce_dimension, pce_multi_index)
    %% Calculating PCE base functions values at quadrature nodes

    % Input:
    %   Q_nodes: Quadrature nodes;
    %   pce_dimension: PCE dimension;
    %   pce_multi_index: Multiple index of PCE;

    % Output:
    %   PCE_base_nodal_value: PCE base functions values at quadrature nodes;

    %%
    PCE_base_nodal_value = zeros(size(Q_nodes, 1), pce_dimension);
    PCE_base_handle = cell(pce_dimension, 1);

    for ii = 1:pce_dimension
        PCE_base_handle{ii} = @(xi_all) pce_base_fun(pce_multi_index(ii, :), xi_all);
    end

    ppm = ParforProgressbar(pce_dimension);

    %%
    parfor ii = 1:pce_dimension
        [PCE_base_nodal_value(:, ii), ~] = f_integration_nodal(PCE_base_handle{ii}, Q_nodes);
        ppm.increment();
    end

    delete(ppm);

    %%
    % save('PCE_base_nodal_value_SG.mat', 'PCE_base_nodal_value');

    %%
    fprintf('\nCalculation for PCE function values at quadrature nodes completed.\n');

end % End of PCE_base_nodal()
