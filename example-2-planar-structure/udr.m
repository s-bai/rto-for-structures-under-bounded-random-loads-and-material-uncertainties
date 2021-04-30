function [c_mean_UDR, Dc_mean_UDR, c_STD_UDR, Dc_STD_UDR, c_mean_MC, Dc_mean_MC, c_STD_MC, Dc_STD_MC] = ...
        udr(mesh_info, rand_numbr_f, eole_order_e, eole_order_t, GH_nodes, GH_weights, u_handle, n_monte_carlo)
    %% Calculating statistical moments using Univariate dimension reduction method
    % Written by Song Bai
    % Supervised by Prof. Zhan Kang

    % Refere to "Univariate Dimension Reduction.docx" for details.

    % Input:
    %       mesh_info.nelx: Num. of elements in direction x;
    %       mesh_info.nely: Num. of elements in direction y;

    %       rand_numbr_f:   Num. of random variables for loads;
    %       eole_order_e:   Num. of random variables for Young's random field;
    %       eole_order_t:   Num. of random variables for thickness random field;

    %       n_tp_nodes_1d: Num. of 1-dimensional base nodes for tensor product grid;
    %       u_handle: Structural response function handel;
    %       n_monte_carlo: Number of Monte Carlo simulation samples;

    % Output:
    %       mean_response_UDR:    Mean value of response calculated by UDR;
    %       mean_sensitivity_UDR: Mean value of response sensitivity calculated by UDR;
    %       std_response_UDR:     STD of response calculated by UDR;
    %       std_sensitivity_UDR:  STD of response sensitivity calculated by UDR;
    %       mean_response_MC:     Mean value of response calculated by Monte Carlo method;
    %       mean_sensitivity_MC:  Mean value of response sensitivity calculated by
    %       std_response_MC:      STD of response calculated by Monte Carlo method;
    %       std_sensitivity_MC:   STD of response sensitivity calculated by Monte Carlo method;

    nelx = mesh_info.nelx;
    nely = mesh_info.nely;

    % Total random variables number
    n_rand_var = rand_numbr_f + eole_order_e + eole_order_t;

    % % Generating one dimensional Gauss-Hermite quadrature nodes
    % [GH_nodes, GH_weights] = gauss_hermite_node(n_tp_nodes_1d);

    %% Quadrature by UDR method
    % Generating quadrature nodes
    UDR_nodes = kron(eye(n_rand_var), GH_nodes);
    % Variable transformation for representation of standard Gaussian random variable
    UDR_nodes = sqrt(2) * UDR_nodes;
    UDR_weights = repmat(GH_weights, n_rand_var, 1);

    % Num. of Gauss-Hermite quadrature nodes
    n_quad_nodes = size(GH_weights, 1);

    % Temp variable storing response value
    response_gi = zeros(n_quad_nodes, 1);

    % Temp variable storing square of response value
    response_gi_square = zeros(n_quad_nodes, 1);

    % Temp variable storing sensitivity value
    response_sensitivity_Dgi = zeros(nely, nelx, n_quad_nodes);

    % Temp variable storing g_i*Dg_i/Drho_e
    temp_gi_Multpy_Dgi = zeros(nely, nelx, n_quad_nodes);

    % Temp variable storing E[gi]*E[Dgi]
    temp_Egi_Multpy_EDgi = zeros(nely, nelx, n_rand_var);

    %%
    xi_all = UDR_nodes;

    %%
    parfor ii = 1:size(UDR_nodes, 1)% The # of columns equals the # of variables

        [~, c, Dc] = feval(u_handle, xi_all(ii, :));

        response_gi(ii, 1) = ...
            1 / (sqrt(pi)) * ...
            c * ...
            UDR_weights(ii, 1);

        response_gi_square(ii, 1) = ...
            1 / (sqrt(pi)) * ...
            c^2 * ...
            UDR_weights(ii, 1);

        response_sensitivity_Dgi(:, :, ii) = ...
            1 / (sqrt(pi)) * ...
            Dc * ...
            UDR_weights(ii, 1);

        temp_gi_Multpy_Dgi(:, :, ii) = ...
            1 / (sqrt(pi)) * ...
            2 * c * Dc * ...
            UDR_weights(ii, 1);
    end

    %%
    for ii = 1:n_rand_var
        E_gi = sum(response_gi((ii - 1) * n_quad_nodes + 1:ii * n_quad_nodes, 1), 1);
        E_Dgi = sum(response_sensitivity_Dgi(:, :, (ii - 1) * n_quad_nodes + 1:ii * n_quad_nodes), 3);

        temp_Egi_Multpy_EDgi(:, :, ii) = 2 * E_gi * E_Dgi;
    end

    xi_all_mean = 0 * UDR_nodes(1, :);
    [~, c_mean, Dc_mean] = feval(u_handle, xi_all_mean);

    %% Calculating mean value of structural response
    c_mean_UDR = sum(response_gi) - (n_rand_var - 1) * c_mean;

    %% Calculating sensitivity of structural response (compliance)
    Dc_mean_UDR = sum(response_sensitivity_Dgi, 3) - (n_rand_var - 1) * Dc_mean;
    %     mean_sensitivity_UDR = sum(response_sensitivity_Dgi, 3);

    %% Calculating STD of structural response (compliance)
    response_square_mean_temp = 0;

    for jj = 1:n_rand_var
        response_square_mean_temp = ...
            response_square_mean_temp + ...
            sum(response_gi(1 + (jj - 1) * size(GH_weights, 1):jj * size(GH_weights, 1), 1))^2;
    end

    c_STD_UDR = sqrt(sum(response_gi_square) - response_square_mean_temp);

    %% Calculating STD sensitivity of structural response (compliance)
    Dc_STD_UDR = 1/2 * (sum(response_gi_square) - response_square_mean_temp)^(-1/2) * ...
        (sum(temp_gi_Multpy_Dgi, 3) - sum(temp_Egi_Multpy_EDgi, 3));

    %% Quadrature by Monte Carlo method
    if n_monte_carlo > 0
        %% Initialization
        response_MC_temp = zeros(n_monte_carlo, 1);
        sensitivity_MC_temp = zeros(nely, nelx, n_monte_carlo);

        c_times_Dc = zeros(nely, nelx, n_monte_carlo);

        %% Generating sample nodes (Standard Gaussian distribution random numbers)
        MC_nodes = randn(n_monte_carlo, n_rand_var);

        %% Calculating mean and STD by Monte Carlo method

        parfor ii = 1:n_monte_carlo
            xi_all = MC_nodes(ii, :);

            [~, c, Dc] = feval(u_handle, xi_all);
            response_MC_temp(ii, 1) = c;
            sensitivity_MC_temp(:, :, ii) = Dc;
            c_times_Dc(:, :, ii) = c * Dc;
        end

        c_mean_MC = mean(response_MC_temp, 1);
        Dc_mean_MC = mean(sensitivity_MC_temp, 3);

        %     std_response_MC = std(response_MC_temp, 0, 1);
        c_STD_MC = std(response_MC_temp);
        Dc_STD_MC = 1 / c_STD_MC * (mean(c_times_Dc, 3) - c_mean_MC * Dc_mean_MC);

    end

end % End of udr()
