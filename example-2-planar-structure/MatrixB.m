function MatrixBvalue = MatrixB(xi_eta)
    %% Matrix B: The strain matrix
    % Input:
    %   xi_eta: Local coordinates

    % Output:
    %   MatrixBvalue: Matrix B at the given points

    % Note:
    %   Here, the normalized Cartesian element is adopted. Hence, the inverse
    %   of Jacobi matrix is
    %                        [2, 0;
    %                         0, 2]
    %   for each element.

    MatrixB1 = 1/2 * [
                (-1 + xi_eta(2)), 0;
                0, (-1 + xi_eta(1));
                (-1 + xi_eta(1)), (-1 + xi_eta(2))];

    MatrixB2 = 1/2 * [
                (1 - xi_eta(2)), 0;
                0, (-1 - xi_eta(1));
                (-1 - xi_eta(1)), (1 - xi_eta(2))];

    MatrixB3 = 1/2 * [
                (1 + xi_eta(2)), 0;
                0, (1 + xi_eta(1));
                (1 + xi_eta(1)), (1 + xi_eta(2))];

    MatrixB4 = 1/2 * [
                (-1 - xi_eta(2)), 0;
                0, (1 - xi_eta(1));
                (1 - xi_eta(1)), (-1 - xi_eta(2))];

    MatrixBvalue = [MatrixB1, MatrixB2, MatrixB3, MatrixB4];

end % End of MatrixB
