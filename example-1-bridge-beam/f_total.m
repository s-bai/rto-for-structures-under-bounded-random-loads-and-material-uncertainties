function [uniform_rand_var, f_random, f_sparse] = f_total(xi_f, f_components, w, nelx, nely)
%% Building total sparse loads vector

% Note: All vectors are stored column-wisely

% Input:
%   xi_f: Underlying Gaussian random variables for defining loads;
%       The odd elements store the random variables for angles;
%       The even elements store the random variables for diameters;

%   w:
%       w(ii).w_matrix: Characteristic matrix of ellipsoid convex model;
%       w(ii).Q: Eigenvectors matrix;
%       w(ii).SIGMA: Eigenvalues diagonal matrix;
%       w(ii).T:
%       w(ii).T_inverse: Inverse of matrix T;

%   f_components:
%       f_components.n: Number of loads;
%       f_components.DOF: DOFs of loads;
%       f_components.nominal_value: Nominal values of loads;

% Output:
%   f_sparse: Sparse vector of loads

%%
n = f_components.n; % Number of loads
DOFs = f_components.DOF; % DOFs of loads
nominal_values = f_components.nominal_value; % Nominal loads values



%% Calculating random variables in normalized psi space
uniform_rand_var = zeros(n, 1);
psi_xy = zeros(n, 1);
mean_mu = 0;
sd_sigma = 1;
cdf_gaussian = cdf('Normal', xi_f, mean_mu, sd_sigma);


%% Calculating random variables in normalized box model
% for ii = 1:n/2
%     uniform_rand_var(2*ii - 1, 1) = unifinv(cdf_gaussian(2*ii - 1), -1.118, 1.118);
%     uniform_rand_var(2*ii, 1) = unifinv(cdf_gaussian(2*ii), 1, 2);
% 
% 
% end
% f_random = zeros(n, 1);
% for ii = 1:n/2
%     f_random(2*ii - 1:2*ii, 1) = [uniform_rand_var(2*ii - 1, 1), uniform_rand_var(2*ii, 1)];
% 
% 
% end


%%
if n == 0
    f_random = f_components.nominal_value;
    
else
    %%
    for ii = 1:n/2
        uniform_rand_var(2*ii - 1, 1) = unifinv(cdf_gaussian(2*ii - 1), -pi/2, pi/2);
        uniform_rand_var(2*ii, 1) = unifinv(cdf_gaussian(2*ii), -1, 1);
        uniform_rand_var(2*ii, 1) = sign(uniform_rand_var(2*ii, 1))*sqrt(abs(uniform_rand_var(2*ii, 1)));
        % uniform_rand_var(2*ii, 1) = sqrt(uniform_rand_var(2*ii, 1));
        
        % Calculating coordinates in psi space
        psi_xy(2*ii - 1, 1) = uniform_rand_var(2*ii, 1)*cos(uniform_rand_var(2*ii - 1, 1));
        psi_xy(2*ii, 1) = uniform_rand_var(2*ii, 1)*sin(uniform_rand_var(2*ii - 1, 1));
    end
    
    %% Transforming loads in psi space to physical space
    f_random = zeros(n/2, 1);
    for ii = 1:n/2
        f_random(2*ii - 1:2*ii, 1) = ...
            w(ii).Q_physical*w(ii).T_interim_inverse*psi_xy(2*ii - 1:2*ii, 1) + nominal_values(2*ii - 1:2*ii, 1);
%         if f_random(2*ii - 1, 1) ~= 0
%             fprintf('\nfx = %.2f\n', f_random(2*ii - 1, 1));
%         end
    end
    
end



%% Building total sparse loads
f_sparse = sparse(DOFs, ones(size(f_components.nominal_value)), f_random, 2*(nelx + 1)*(nely + 1), 1);


end % End of f_total()



