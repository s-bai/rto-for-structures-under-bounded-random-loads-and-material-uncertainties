function [zeta_value, delta_value] = random_field_fun(xi, phi, mesh_info, uniform_dist_boundary)
%% Calculating the random field

% Input:
%   xi: The uncorrelated standard Gaussian random variables;
%   phi: Mode functions;
%   mesh_info:
%       mesh_info.nelx: Num. of elements in direction x;
%       mesh_info.nely: Num. of elements in direction y;
%       mesh_info.ik: The row index vectors
%       mesh_info.jk: The column index vectors
%       mesh_info.element_node_num: Element node numbers;
%   uniform_dist_boundary: Upper and lower boundaries of the uniform distribution.

% Output:
%   zeta: The underlying Gaussian random field (In matrix form);
%   delta: The perturbation random field for physical quantities (In matrix form).
%
% Note:
%   Sample points of the random field are sampled at nodes.

%% Num. of elements in directions x and y
nelx = mesh_info.nelx;
nely = mesh_info.nely;

%% Calculating the underlying Gaussian random field
zeta_value = reshape(sum(xi.'.*phi), [nely + 1, nelx + 1]);

%% Calculating the perturbation random field for physical quantities
mean_mu = 0;
sd_sigma = 1;
cdf_value = cdf('Normal', zeta_value, mean_mu, sd_sigma);
delta_value = unifinv(cdf_value, uniform_dist_boundary(1), uniform_dist_boundary(2));
