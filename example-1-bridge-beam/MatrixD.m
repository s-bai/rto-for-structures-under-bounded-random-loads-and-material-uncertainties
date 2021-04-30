function [matrix_D, d_matrix_D] = MatrixD(local_coordinate, E_Ke, t_Ke, nu_Ke, rho, penalty)
%% Matrix D: THe elasticity matrix

% Input:
%   local_coordinate: Local coordinates of Gaussian points
%   E_Ke: Elemental Young's modulus;
%   t_Ke: Elemental thickness;
%   nu_Ke: Elemental Poisson's ratio;
%   rho: Element density;
%   penalty: Penalty factor;

% Output:
%   matrix_D: Elasticity matrix value at the given point
%   d_matrix_D: Derivative of matrix D with respect to element density rho

STIFFNESS_MIN = 1e-9;

nodal_element_stiff_factor = E_Ke.*t_Ke;
nodal_element_stiff_factor_penalized = ...
    STIFFNESS_MIN + (nodal_element_stiff_factor - STIFFNESS_MIN).*rho^penalty;

element_stiff_factor = sum(ShapeFun(local_coordinate).*nodal_element_stiff_factor_penalized);

matrix_D = element_stiff_factor/(1 - nu_Ke^2)*...
    [
    1,     nu_Ke, 0;
    nu_Ke, 1,     0;
    0,     0,     (1 - nu_Ke)/2];

if nargout > 1
    d_nodal_element_stiff_factor_penalized = ...
        penalty*(nodal_element_stiff_factor - STIFFNESS_MIN).*rho^(penalty - 1);
    
    d_element_stiff_factor = sum(ShapeFun(local_coordinate).*d_nodal_element_stiff_factor_penalized);
    
    d_matrix_D = d_element_stiff_factor/(1 - nu_Ke^2)*...
        [
        1,     nu_Ke, 0;
        nu_Ke, 1,     0;
        0,     0,     (1 - nu_Ke)/2];
end

end % End of MatrixD