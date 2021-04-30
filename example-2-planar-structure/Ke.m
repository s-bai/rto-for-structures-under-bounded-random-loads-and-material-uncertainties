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
    G = 0.577350269189626;
    GaussQuadPts = transpose([-G, -G; G, -G; G, G; -G, G]);
    % Weights of Gaussian quadrature points
    GaussQuadWeights = ones(2, 4);
    % Determinant of Jacobi matrix
    JACOBI_DET = 1/4;

    Ke_value = zeros(8, 8);

    for ii = 1:4
        Ke_value = Ke_value + ...
            MatrixB(GaussQuadPts(:, ii))' * ...
            MatrixD(GaussQuadPts(:, ii), e_Ke, t_Ke, nu_Ke, rho, penalty) * ...
            MatrixB(GaussQuadPts(:, ii)) * ...
            JACOBI_DET * ...
            prod(GaussQuadWeights(:, ii));
    end

    if nargout > 1

        d_Ke_value = zeros(8, 8);

        for ii = 1:4
            [~, d_matrix_D] = MatrixD(GaussQuadPts(:, ii), e_Ke, t_Ke, nu_Ke, rho, penalty);

            d_Ke_value = d_Ke_value + ...
                MatrixB(GaussQuadPts(:, ii))' * ...
                d_matrix_D * ...
                MatrixB(GaussQuadPts(:, ii)) * ...
                JACOBI_DET * ...
                prod(GaussQuadWeights(:, ii));
        end

    end

end % End of Ke()
