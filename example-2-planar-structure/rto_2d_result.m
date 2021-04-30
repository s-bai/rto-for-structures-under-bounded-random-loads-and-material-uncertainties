function rto_2d_result(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid, n_MC_samples)
    %% Plotting iteration history 
    % Input:
    %	See input explaination in the main function.

    % Output:
    %	Iteration histories.

    %% Load result file
    cmd_name = 'rto_2d';
    cmd_parameters = [nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid];

    iteration_file_name = cmd_history(cmd_name, cmd_parameters);

    %% Load RDO results
    % 	'rho_history'
    % 	'objective_iter_history'
    %   'constraint_iter_history'
    %   'iter'
    % 	'n_ellipse'
    rdo_iter_history = load(iteration_file_name);

    % try
    %     n_ellipse = rdo_iter_history.n_ellipse;
    % catch
    %     prompt = '# of ellipse model is missing, manually input:';
    %     n_ellipse = input(prompt);
    % end

    try
        max_iter = rdo_iter_history.iter;
    catch
        temp = ...
            find(rdo_iter_history.worst_case_compliance_history(:, 1));
        max_iter = temp(end);
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Boundary conditions
    % boundary_condition = rdo_iter_history.boundary_condition;
    %
    % % Problem paramaters
    % nu = rdo_iter_history.problem_parameters.nu;
    %
    % e = rdo_iter_history.problem_parameters.e*ones(nely + 1, nelx + 1);
    % err_magnitude_e = rdo_iter_history.problem_parameters.err_magnitude_e;
    % correlation_length_e = rdo_iter_history.problem_parameters.correlation_length_e;
    %
    % t = rdo_iter_history.problem_parameters.t*ones(nely + 1, nelx + 1);
    % err_magnitude_t = rdo_iter_history.problem_parameters.err_magnitude_t;
    % correlation_length_t = rdo_iter_history.problem_parameters.correlation_length_t;
    %
    % eole_order_e = rdo_iter_history.problem_parameters.eole_order_e;
    % eole_order_t = rdo_iter_history.problem_parameters.eole_order_t;
    %
    % ellips_char_parameter = rdo_iter_history.ellips_char_parameter;
    %
    % uniform_dist_boundary_e = rdo_iter_history.uniform_dist_boundary_e;
    % uniform_dist_boundary_t = rdo_iter_history.uniform_dist_boundary_t;
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %% Load deterministic results
    % % Note that the file name of deterministic results is input manually
    % cmd_name = 'det_rto_against_uncertainty';
    % cmd_parameters = [nelx,nely,volfrac,penalty,rmin];
    % iteration_file_name = cmd_history(cmd_name, cmd_parameters);
    % det_iter_history = load(iteration_file_name);

    %% Plotting iteration histories of mean and std
    figure('Name', 'Iteration history of mean and std', 'NumberTitle', 'off');

    n_iter_steps = linspace(1, max_iter, max_iter);

    % Plotting iteration history of mean value
    yyaxis left;
    plot(n_iter_steps, rdo_iter_history.objective_iter_history(1:max_iter, 3));
    ylim([0, 1.5 * max(rdo_iter_history.objective_iter_history(1:max_iter, 3))]);

    % Plotting iteration history of STD
    yyaxis right;
    plot(n_iter_steps, rdo_iter_history.objective_iter_history(1:max_iter, 4));
    ylim([0, 1.5 * max(rdo_iter_history.objective_iter_history(1:max_iter, 4))]);

    %% Plotting optimal topology layout
    figure('Name', 'Robust optimal topology layout', 'NumberTitle', 'off');

    optimal_topology = rdo_iter_history.rho_history(:, :, max_iter + 1);
    % optimal_topology = rdo_iter_history.rho_history(:, :, 10);

    colormap(gray);
    imagesc(1 - optimal_topology);
    caxis([0 1]);
    axis equal;
    axis off;
    drawnow;

    %% Comparing mean and std obtained by PCE and UDR and MC

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Boundary conditions
    boundary_condition = rdo_iter_history.boundary_condition;

    % Problem paramaters
    nu = rdo_iter_history.problem_parameters.nu;

    e = rdo_iter_history.problem_parameters.e * ones(nely + 1, nelx + 1);
    err_magnitude_e = rdo_iter_history.problem_parameters.err_magnitude_e;
    correlation_length_e = rdo_iter_history.problem_parameters.correlation_length_e;

    t = rdo_iter_history.problem_parameters.t * ones(nely + 1, nelx + 1);
    err_magnitude_t = rdo_iter_history.problem_parameters.err_magnitude_t;
    correlation_length_t = rdo_iter_history.problem_parameters.correlation_length_t;

    eole_order_e = rdo_iter_history.problem_parameters.eole_order_e;
    eole_order_t = rdo_iter_history.problem_parameters.eole_order_t;

    % ellips_char_parameter = rdo_iter_history.ellips_char_parameter;

    uniform_dist_boundary_e = rdo_iter_history.uniform_dist_boundary_e;
    uniform_dist_boundary_t = rdo_iter_history.uniform_dist_boundary_t;

    %% f
    f_components = rdo_iter_history.f_components;

    %%
    rand_numbr_f = rdo_iter_history.rand_numbr_f;

    %% w
    w = rdo_iter_history.w;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Plotting PDF and CDF of compliance of optimal topology layout

    %% FEA setup
    % Meshing
    mesh_info = Meshing(nelx, nely);

    % Calculating mode functions of the Young's modulus random field
    if eole_order_e > 0
        [phi_e] = eole_mode_fun(mesh_info, correlation_length_e, eole_order_e);
    else
        phi_e = 0;
    end

    % Calculating mode functions of the thickness random field
    if eole_order_t > 0
        [phi_t] = eole_mode_fun(mesh_info, correlation_length_t, eole_order_t);
    else
        phi_t = 0;
    end

    % Nominal loads
    %   The first column stores the DOFs of the loads
    %   The second column stores the magnitudes
    % Loading at right mid-point
    % problem_parameters.f_nominal = [
    %     2*nelx*(nely + 1) + nely + 1,  1;
    %     2*nelx*(nely + 1) + nely + 2,  0];

    %% Note %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % The nominal loads are saved in the data file as problem_parameters.f_DOF and
    % % f_nominal_val.
    % f_DOF = rdo_iter_history.problem_parameters.f_DOF;
    % f_nominal_val = rdo_iter_history.problem_parameters.f_nominal_val;

    %% Load deterministic results
    % Note that the file name of deterministic results is input manually
    % cmd_name = 'det_rto_2d';
    % cmd_parameters = [nelx, nely, volfrac, penalty, rmin];
    % iteration_file_name = cmd_history(cmd_name, cmd_parameters);
    % det_iter_history = load(iteration_file_name);

    % Set topology layout
    rho_rdo = optimal_topology;
    % rho_det = det_iter_history.x_history(:, :, det_iter_history.iter + 1);

    % figure('Name', 'Deterministic optimal topology layout', 'NumberTitle', 'off');
    %
    % colormap(gray);
    % imagesc(1 - rho_det);
    % caxis([0 1]);
    % axis equal;
    % axis off;
    % drawnow;

    prompt_text = 'Proceed to calculate STD of total volume of deterministic design (y/n)?';
    switch_str = input(prompt_text, 's');

    is_det_vol_std = 0;

    if strcmp(switch_str, 'y')
        %% Calculating mean and STD of total volume of deterministic design
        is_det_vol_std = 1;
        prompt = 'Input pce_order_t = ';
        pce_order_t = input(prompt);

        [volume_pce_dimension, volume_pce_multi_index] = generate_pce_multi_index(eole_order_t, pce_order_t);

        pce_gamma_volume = zeros(volume_pce_dimension, 1);

        for ii = 1:volume_pce_dimension
            pce_gamma_volume(ii, 1) = pce_normal_factor(volume_pce_multi_index(ii, :));
        end

        v_hat = zeros(volume_pce_dimension, 1);
        det_volume_STD = zeros(det_iter_history.iter + 1, 1);

        det_ii = 1;

        fprintf('\nVolume verification of deterministic design progress: ');

        while det_ii <= (det_iter_history.iter + 1)

            det_rho_temp = det_iter_history.x_history(:, :, det_ii);

            % Function handle for volume response
            v_handle = @(xi_t) total_volume(mesh_info, t, det_rho_temp, xi_t, phi_t, err_magnitude_t, uniform_dist_boundary_t);

            % Function handle for PCE base function
            for jj = 1:volume_pce_dimension
                v_pce_base_handle = @(xi_t) ...
                    pce_base_fun(volume_pce_multi_index(jj, :), 0, 0, eole_order_t, xi_t);

                [response_integration_sg, ~, ~] = ...
                    gh_integration(mesh_info, 0, 0, eole_order_t, ...
                    k_level, n_tp_grid, v_handle, v_pce_base_handle);

                v_hat(jj) = (1 / pce_gamma_volume(jj)) * response_integration_sg;

            end

            det_volume_STD(det_ii, 1) = sqrt(sum(pce_gamma_volume(2:end) .* (v_hat(2:end).^2)));

            if det_ii > 1

                for temp = 1:6
                    fprintf('\b');
                end

            end

            fprintf('%5.1f%%', det_ii / (det_iter_history.iter + 1) * 100);
            det_ii = det_ii + 1;

        end

        fprintf('\n')

    end

    if is_det_vol_std == 1
        %% Saving det_volume_STD
        save('det_volume_STD.mat', 'det_volume_STD');
    end

    %% Plotting compliance iteration history and error bar of volume constraint of RDO
    prompt = 'Proceed to plot volume contraint history obtained by PCE (1) or UDR (2)?';
    is_pce_udr = input(prompt);

    switch is_pce_udr
        case 1
            is_pce_udr = [1, 2];
        case 2
            is_pce_udr = [3, 4];
    end

    if eole_order_t > 0
        %     load rdo results
        % load det results
        % load rdo iter history of objective (mean and std of compliance)
        % 						mean and std of volume constraint
        % load det iter history of objective (compliance)
        % 						deterministic volume fraction value and the corresponding std
        % plotting objective iter history	and error bar (STD) of rdo design in the same figure using "hold on" command
        %
        % plotting the counterpart of the above results of deterministic design

        figure('Name', 'Iteration history of compliance and error bar of volume fraction', 'NumberTitle', 'off');

        n_iter_steps = linspace(1, max_iter, max_iter);

        % Plotting iteration history of mean value
        yyaxis left;

        plot(n_iter_steps, rdo_iter_history.objective_iter_history(1:max_iter, is_pce_udr(1)));
        ylim([0, 1.5 * max(rdo_iter_history.objective_iter_history(1:max_iter, is_pce_udr(1)))]);

        % Plotting error bar of volume fraction
        volume_mean = rdo_iter_history.constraint_iter_history(:, 1);
        volume_STD = rdo_iter_history.constraint_iter_history(:, 2);
        volume_fraction_mean = volume_mean / (nelx * nely);
        volume_fraction_STD = volume_STD / (nelx * nely);

        yyaxis right;
        plot(n_iter_steps, volume_fraction_mean(1:max_iter, 1));
        errorbar(n_iter_steps, volume_fraction_mean(1:max_iter, 1), volume_fraction_STD(1:max_iter, 1));
        ylim([0, 1]);

        %% Drawing boundedline
        figure;
        [line_handle, patch_handle] = ...
            boundedline(n_iter_steps, volume_fraction_mean(1:max_iter, 1), volume_fraction_STD(1:max_iter, 1), '-b');
        outlinebounds(line_handle, patch_handle);
        title('Mean and STD history of volume fraction');
        axis tight;
        ylim([0, 1]);

    end

    %% Plotting compliance iteration history and error bar of volume constraint of DET
    prompt_text = 'Proceed to draw volume fraction error bar of deterministic design? (y/n)';
    switch_str = input(prompt_text, 's');

    if strcmp(switch_str, 'y')

        try
            det_vol_fraction_history = zeros((det_iter_history.iter + 1), 1);

            for ii = 1:(det_iter_history.iter + 1)
                det_vol_fraction_history(ii, 1) = ...
                    sum(sum(det_iter_history.x_history(:, :, ii))) / (nelx * nely);
            end

            det_volume_STD = load('det_volume_STD.mat');
            det_volume_STD_history = zeros((det_iter_history.iter + 1), 1);

            for ii = 1:(det_iter_history.iter + 1)
                det_volume_STD_history(ii, 1) = ...
                    det_volume_STD.det_volume_STD(ii) / (nelx * nely);
            end

            %% Drawing boundedline
            figure;
            det_total_iter = det_iter_history.iter + 1;
            [line_handle, patch_handle] = ...
                boundedline(1:det_total_iter, det_vol_fraction_history(1:det_total_iter, 1), det_volume_STD_history(1:det_total_iter, 1), '-b');
            outlinebounds(line_handle, patch_handle);
            title('Mean and STD history of volume fraction');
            axis tight;
            ylim([0, 1]);

        catch
            fprintf('File(s) missing.\n');
        end

    end

    %% Generating function handle for FEA (new)
    u_handle = @(xi_all, rho) FEA(mesh_info, e, eole_order_e, t, eole_order_t, nu, boundary_condition, ...
        f_components, rand_numbr_f, rho, penalty, w, ...
        phi_e, err_magnitude_e, uniform_dist_boundary_e, ...
        phi_t, err_magnitude_t, uniform_dist_boundary_t, ...
        xi_all);

    %% Monte Carlo simulation
    % prompt_text = 'Proceed to Monte Carlo simulation? (y/n)';
    % switch_str = input(prompt_text, 's');
    %
    % if strcmp(switch_str,'y')
    %% Generating Monte Carlo simulation sample nodes
    n_rand_var = rand_numbr_f + eole_order_e + eole_order_t;

    if n_MC_samples > 0
        n_monte_carlo_samples = n_MC_samples;
    else
        n_monte_carlo_samples = 100;
    end

    MC_nodes = randn(n_monte_carlo_samples, n_rand_var);

    % %% ____________UDR_________________
    % [mean_response_UDR, ~, ...
    %     std_response_UDR, ~, ...
    %     mean_response_MC, ~, ...
    %     std_response_MC, ~] = ...
    %     udr(mesh_info, n_rand_var, 0, 0, 5, u_handle, n_MC_samples);

    %% Monte Carlo simulation
    c_samples_rdo = zeros(n_monte_carlo_samples, 1);
    c_samples_det = zeros(n_monte_carlo_samples, 1);

    data_queue_temp = parallel.pool.DataQueue;
    echo_message = waitbar(0, 'Monte Carlo simulation in progress');
    afterEach(data_queue_temp, @nUpdateWaitbar);
    progress_percentage = 1;

    %% Calculating mean and STD by Monte Carlo method (new)
    parfor ii = 1:n_monte_carlo_samples
        xi_all = MC_nodes(ii, :);

        [~, c_rdo, ~] = feval(u_handle, xi_all, rho_rdo);
        c_samples_rdo(ii, 1) = c_rdo;

        %         [~, c_det, ~] = feval(u_handle, xi_all, rho_det);
        %         c_samples_det(ii, 1) = c_det;

        send(data_queue_temp, ii);
    end

    function nUpdateWaitbar(~)
        waitbar(progress_percentage / n_monte_carlo_samples, echo_message);
        progress_percentage = progress_percentage + 1;
    end

    [c_pdf_rdo, c_xi_rdo] = ksdensity(c_samples_rdo, 'Support', 'positive');
    %     [c_pdf_det, c_xi_det]=ksdensity(c_samples_det, 'Support','positive');

    mean_rdo = mean(c_samples_rdo, 1);
    std_rdo = std(c_samples_rdo);

    %     mean_det = mean(c_samples_det, 1);
    %     std_det = std(c_samples_det);

    fprintf('Mean by RDO = %.4f\n', mean_rdo);
    fprintf('STD by RDO = %.4f\n', std_rdo);
    %     fprintf('Mean by DET = %.4f\n', mean_det);
    %     fprintf('STD by DET = %.4f\n', std_det);

    figure;
    plot(c_xi_rdo, c_pdf_rdo, 'LineWidth', 2, 'Color', 'b');
    hold on
    %     plot(c_xi_det, c_pdf_det, 'LineWidth', 2, 'Color', 'r');

    %     legend('Robust design', 'Deterministic design');
    legend('Robust design');
    hold off

    % % objective_iter_history: 1st column stores Mean value iteration history
    % % objective_iter_history: 2nd column stores Standard deviation iteration history
    % objective_iter_history = zeros(max_iter_step, 6);
    %
    % % constraint_iter_history: Non-zero when randomness in thickness is considered
    % constraint_iter_history = zeros(max_iter_step, 2);

    % errorbar

    % end

end % End 
