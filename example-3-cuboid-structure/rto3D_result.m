function rto3D_result(nelx, nely, nelz, n_MC_sample)
    %% Monte Carlo simulation for deterministic design and robust design
    % Input:
    %   n_MC_sample: Number of Monte Carlo simulation samples;

    % Output:
    %   Saving 'c_det' and 'c_rto' to rto_loads3D_MC_results.mat

    %% Load deterministic design results
    det_result = load('det3D.mat');
    rho_det = det_result.rho_history(:, :, :, det_result.iter);

    %% Load robust design results
    rto_result = load('rto3D.mat');
    rho_rto = rto_result.rho_history(:, :, :, rto_result.iter);

    %% Generating Monte Carlo simulation samples
    MC_samples = randn(n_MC_sample, 3);

    %% Initialization for FEA
    % Loading point coordinates: Cube beam with concentrate loads +Fx at right center
    [il, jl, kl] = meshgrid(nelx, nely / 2, nelz / 2);

    % BCs coordinates: Cantilever beam fixed at left facet
    [iif, jf, kf] = meshgrid(0, 0:nely, 0:nelz);

    % Material properties
    E0 = 1; % Young's modulus of solid material
    Emin = 1e-9; % Young's modulus of void-like material
    nu = 0.3; % Poisson's ratio

    % Parameters for uncertain loads
    w = [3, 0, 0; 0, 3, 0; 0, 0, 1];
    f0 = [1; 0; 0];

    % FEA parameters
    [mesh_info, material_properties] = FEA3D_setup(nelx, nely, nelz, il, jl, kl, iif, jf, kf, E0, Emin, nu, w, f0);

    %% Initialization for Monte Carlo simulation
    c_det = zeros(n_MC_sample, 1);
    c_rto = zeros(n_MC_sample, 1);

    FEA_handle = @(rho, xi_all) FEA3D(material_properties, mesh_info, rho, penal, f_components, xi_all);

    %%
    parfor ii = 1:n_MC_sample
        [~, c_det(ii, 1), ~] = feval(FEA_handle, rho_det, MC_samples(ii, :));
        [~, c_rto(ii, 1), ~] = feval(FEA_handle, rho_rto, MC_samples(ii, :));
    end

    %%
    c_mean_det = mean(c_det, 1);
    c_STD_det = std(c_det);

    c_mean_rto = mean(c_rto, 1);
    c_STD_rto = std(c_rto);

    fprintf('\nMean by DET = %.4f, STD by DET = %.4f\n', c_mean_det, c_STD_det);
    fprintf('\nMean by RTO = %.4f, STD by RTO = %.4f\n', c_mean_rto, c_STD_rto);

    %%
    [c_pdf_det, c_xi_det] = ksdensity(c_det, 'Support', 'positive');
    [c_pdf_rto, c_xi_rto] = ksdensity(c_rto, 'Support', 'positive');

    %%
    figure;
    plot(c_xi_rto, c_pdf_rto, 'LineWidth', 2, 'Color', 'b');
    hold on
    plot(c_xi_det, c_pdf_det, 'LineWidth', 2, 'Color', 'r');

    legend('Robust design', 'Deterministic design');
    hold off

    %%
    save('rto_loads3D_MC_results.mat', 'c_det', 'c_rto');

end % End of rto3D_result()
