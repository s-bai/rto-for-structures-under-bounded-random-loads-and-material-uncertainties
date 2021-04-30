function [zeta_value, delta_value] = gaussian_to_uniform(xi_all_e, phi3D)


%%
zeta_value = xi_all_e*phi3D;

%% Calculating the perturbation random field for physical quantities
mean_mu = 0;
std_sigma = 1;

cdf_value = cdf('Normal', zeta_value, mean_mu, std_sigma);
delta_value = unifinv(cdf_value, -1, 1);


end % End of gaussian_to_uniform