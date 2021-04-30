function [Ke_value, d_Ke_value] = Ke(e_Ke, t_Ke, nu_Ke, rho, penalty)
    %% Calculating element stiffness matrix

    % Input:
    %   e_Ke: Element Young's modulus;
    %   nu_Ke: Element Poisson's ratio;
    %   rho: Element density;
    %   t_Ke: Element thickness

    % Output:
    %   Ke_value: Element stiffness matrix
    %   d_Ke_value: Sensitivity of element stiffness matrix with respect to element density rho

    %% Numerical integration parameters setup
    % Gaussian quadrature points
    % G = 0.577350269189626;
    % GaussQuadPts = transpose([-G, -G; G, -G; G, G; -G, G]);
    gauss_pt = [-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626; ...
                -0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626];

    % Weights of Gaussian quadrature points
    GaussQuadWeights = ones(2, 4);
    % Determinant of Jacobi matrix
    JACOBI_DET = 1/4;

    Ke_value = zeros(8, 8);

    for ii = 1:4
        Ke_value = Ke_value + ...
            MatrixB(gauss_pt(:, ii))' * ...
            MatrixD(gauss_pt(:, ii), e_Ke, t_Ke, nu_Ke, rho, penalty) * ...
            MatrixB(gauss_pt(:, ii)) * ...
            JACOBI_DET * ...
            prod(GaussQuadWeights(:, ii));
    end

    if nargout > 1

        d_Ke_value = zeros(8, 8);

        for ii = 1:4
            [~, d_matrix_D] = MatrixD(gauss_pt(:, ii), e_Ke, t_Ke, nu_Ke, rho, penalty);

            d_Ke_value = d_Ke_value + ...
                MatrixB(gauss_pt(:, ii))' * ...
                d_matrix_D * ...
                MatrixB(gauss_pt(:, ii)) * ...
                JACOBI_DET * ...
                prod(GaussQuadWeights(:, ii));
        end

    end

end % End of Ke()
