function [f_random, f_sparse] = f_total3d(xi_f, f_components)
%% Building total sparse loads vector

% Note: All vectors are stored column-wisely

% Input:
%   xi_f: Underlying Gaussian random variables for defining loads;
%       The odd elements store the random variables for angles;
%       The even elements store the random variables for diameters;

%   w:
%       w(ii).w_matrix: Characteristic matrix of ellipsoid convex model;
%       w(ii).Q_physical: Eigenvectors matrix of loads in physical space;
%       w(ii).SIGMA_physical: Eigenvalues diagonal matrix of loads in physical space;
%       w(ii).T_physical_inverse: Inverse of matrix T for loads in physical space;
%       w(ii).T_interim_inverse: Transformation matrix for normalized model to physical model transformation;

%   f_components:
%       f_components.n: Number of loads;
%       f_components.DOF: DOFs of loads;
%       f_components.nominal_value: Nominal values of loads;

% Output:
%   f_sparse: Sparse vector of loads

%%
% Number of loads
n = f_components.n;
% Characteristic matrix
% w = f_components.w;
% 
% Q = f_components.Q;
% SIGMA = f_components.SIGMA;
% T = f_components.T;
T_original_inverse = f_components.T_original_inverse;
% T_interim_inverse = f_components.T_interim_inverse;
% Q_original = f_components.Q_original;

%

% DOFs of loads
loaddof = f_components.DOF;
% Nominal loads
f0 = f_components.f0;
% Total DOFs
ndof = f_components.total_DOFs;





%% Calculating random variables in normalized psi space
uniform_rand_var = zeros(3*n, 1);
psi_xy = zeros(3*n, 1);
mean_mu = 0;
sd_sigma = 1;
cdf_gaussian = cdf('Normal', xi_f, mean_mu, sd_sigma);


% %% Calculating random variables in normalized box model
% for ii = 1:n
%     uniform_rand_var(2*ii - 1, 1) = unifinv(cdf_gaussian(2*ii - 1), -1, 1);
%     uniform_rand_var(2*ii, 1) = unifinv(cdf_gaussian(2*ii), -3, -1);
%
%
% end
% f_random = zeros(2*n, 1);
% for ii = 1:n
%     f_random(2*ii - 1:2*ii, 1) = [uniform_rand_var(2*ii - 1, 1), uniform_rand_var(2*ii, 1)];
%
%
% end


%%
if n == 0
    f_random = f0;
    
else
    %%
    for ii = 1:n
        uniform_rand_var(3*ii - 2, 1) = unifinv(cdf_gaussian(3*ii - 2), 0, pi);
        uniform_rand_var(3*ii - 1, 1) = unifinv(cdf_gaussian(3*ii - 1), -1, 1);
        uniform_rand_var(3*ii,     1) = unifinv(cdf_gaussian(3*ii), -1, 1);
        
        uniform_rand_var(3*ii - 1, 1) = acos(uniform_rand_var(3*ii - 1, 1)); 
        uniform_rand_var(3*ii,     1) = nthroot(uniform_rand_var(3*ii, 1), 3);
        
                
        % Calculating coordinates in psi space
        psi_xy(3*ii - 2, 1) = uniform_rand_var(3*ii, 1)*cos(uniform_rand_var(3*ii - 2, 1))*sin(uniform_rand_var(3*ii - 1, 1));
        psi_xy(3*ii - 1, 1) = uniform_rand_var(3*ii, 1)*sin(uniform_rand_var(3*ii - 2, 1))*sin(uniform_rand_var(3*ii - 1, 1));
        psi_xy(3*ii    , 1) = uniform_rand_var(3*ii, 1)*cos(uniform_rand_var(3*ii - 1, 1));
    end
    
    %% Transforming loads in psi space to physical space
    f_random = zeros(3*n, 1);
    for ii = 1:n
%         f_random(3*ii - 2:3*ii, 1) = T_original_inverse*psi_xy(3*ii - 2:3*ii, 1) + f0(3*ii - 2:3*ii, 1);
        f_random(3*ii - 2:3*ii, 1) = f_components.T*psi_xy(3*ii - 2:3*ii, 1) + f0(3*ii - 2:3*ii, 1);
    end
    
end



%% Building total sparse loads
f_sparse = sparse(loaddof, 1, f_random, ndof, 1);


end % End of f_total()



