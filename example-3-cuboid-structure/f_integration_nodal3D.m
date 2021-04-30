function [f_nodal_value, Df_nodal_value] = f_integration_nodal3D(f_handle, nodes, nelx, nely, nelz)
    %% Function values at numerical integration nodes
    % Input:
    %   f_handle: Handle of function f;
    %   nodes: Quadrature nodes;
    %   nelx, nely, nelz: Division number in x, y, z directions;

    % Output:
    %   f_nodal_value: Values of function f at integration nodes;
    %   Df_nodal_value: Values of function Df at integration nodes;

    %%
    switch nargin
        case 2
            f_nodal_value = zeros(size(nodes, 1), 1);
            Df_nodal_value = 0;

            parfor ii = 1:size(nodes, 1)
                [~, f_nodal_value(ii, 1), ~] = feval(f_handle, nodes(ii, :));
            end

        case 5
            f_nodal_value = zeros(size(nodes, 1), 1);
            Df_nodal_value = zeros(nely, nelx, nelz, size(nodes, 1));

            parfor ii = 1:size(nodes, 1)
                [~, f_nodal_value(ii, 1), Df_nodal_value(:, :, :, ii)] = feval(f_handle, nodes(ii, :));
            end

    end

end % End of f_integration()
