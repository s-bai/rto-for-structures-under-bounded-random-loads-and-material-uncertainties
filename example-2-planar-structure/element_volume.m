function [volume, d_volume] = element_volume(t_Ke, rho)
    %% Calculating total volume

    % Input:
    %   t_Ke: Element thickness
    %   rho: Element density;

    % Output:
    %   volume: Total element volume
    %   d_volume: Sensitivity of element volume with respect to element density.

    % Note:
    %   In the present finite element meshing, the Jacobian determinant equals 1/4.

    %% Numerical integration parameters setup
    % Gaussian quadrature points
    G = 0.577350269189626;
    GaussQuadPts = transpose([-G, -G; G, -G; G, G; -G, G]);
    % Weights of Gaussian quadrature points
    GaussQuadWeights = ones(2, 4);
    % Determinant of Jacobi matrix
    JACOBI_DET = 1/4;

    %% Ininitalization
    volume = 0;
    d_volume = 0;

    for ii = 1:4
        [node_volume_temp, d_node_volume_temp] = node_volume(GaussQuadPts(:, ii), t_Ke, rho);

        volume = volume + ...
            node_volume_temp * ...
            JACOBI_DET * ...
            prod(GaussQuadWeights(:, ii));

        d_volume = d_volume + ...
            d_node_volume_temp * ...
            JACOBI_DET * ...
            prod(GaussQuadWeights(:, ii));
    end

end % End of element_volume()

%%  ---- Local subroutines (Calculating volume on Gaussian quadrature nodes) -------
function [node_volume_value, d_node_volume_value] = node_volume(local_coordinate, t_Ke, rho)

    STIFFNESS_MIN = 1e-9;

    node_element_volume_temp = STIFFNESS_MIN + (t_Ke - STIFFNESS_MIN) .* rho;

    node_volume_value = sum(ShapeFun(local_coordinate) .* node_element_volume_temp);

    d_node_element_volume_temp = (t_Ke - STIFFNESS_MIN);

    d_node_volume_value = sum(ShapeFun(local_coordinate) .* d_node_element_volume_temp);

end % End of node_volume()
