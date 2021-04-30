function [delta] = gaussian_uniform(xi, uniform_limits)
    %% Memoryless trasformation from Gaussian random variable to uniform random variable
    %   Input:
    %       xi: The uncorrelated underlying standard Gaussian random variables;
    %       uniform_limits: Boundaries of uniform distribution;

    %   Output:
    %       delta: The uniform distribution random variable;

    %%
    mean_mu = 0;
    STD_sigma = 1;

    cdf_value = cdf('Normal', xi, mean_mu, STD_sigma);
    delta = unifinv(cdf_value, uniform_limits(1), uniform_limits(2));

end % End of gaussian_uniform()
