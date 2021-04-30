function rto_3D(nelx, nely, nelz, volfrac, penal, rmin, n_tp_grid, n_MC_samples, sn)
    %% Robust topology optimization for 3D structures considering uncertainties
    % Written by Song Bai, supervised by Zhan Kang.
    % The topology optimization code for 3D structures is based on the work by Liu and Tovar

    % Input:
    %   nelx, nely, nelz: Element division in x, y and z directions;
    %   volfrac: Volume fraction;
    %   penal: Penalty factor;
    %   rmin: Minimum filter radius;
    %   sn: Serial number for batch submission.

    %% Prepare file for saving data
    cmd_name = 'rto_3D';
    cmd_parameters = [nelx, nely, nelz, volfrac, penal, rmin, sn];
    cmd_output = cmd_history(cmd_name, cmd_parameters);

    %% Load data from Excel file
    data_table = readtable('rto_3D.xlsx');

    err_magnitude_e = data_table.err_magnitude_e(sn);
    E0 = data_table.E0(sn);
    eole_order_e = data_table.eole_order_e(sn);
    correlation_length_factor = data_table.correlation_length_factor(sn);
    k_level = data_table.k_level(sn);
    pce_order = data_table.pce_order(sn);
    beta0 = data_table.beta0(sn);
    w_filename = strcat('w', num2str(sn), '.xlsx');
    w = readmatrix(w_filename);

    %% Switch for PCE and UDR (0 for false, 1 for true)
    switch_PCE = 1;
    switch_UDR = 0;

    %% USER-DEFINED LOOP PARAMETERS
    max_iter = 500; % Maximum number of iterations
    min_change = 0.001; % Termination criterion
    % displayflag = 0;  % Display structure flag

    %% USER-DEFINED MATERIAL PROPERTIES
    Emin = 1e-9; % Young's modulus of void-like material
    nu = 0.3; % Poisson's ratio

    %% PREPARE FINITE ELEMENT ANALYSIS
    nele = nelx * nely * nelz;
    ndof = 3 * (nelx + 1) * (nely + 1) * (nelz + 1);

    %% Loads: Cube with concentrate loads +Fx at right center
    % % Coordinates
    % [il,jl,kl] = meshgrid(nelx, nely/2, nelz/2);
    % % Node IDs
    % loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);
    % % DOFs
    % loaddof = [3*loadnid(:) - 2, 3*loadnid(:) - 1, 3*loadnid(:)];

    %% Loads: 3D wheel
    % Coordinates
    il = nelx; jl = nely / 2; kl = nelz / 2;
    % Node IDs
    loadnid = kl * (nelx + 1) * (nely + 1) + il * (nely + 1) + (nely + 1 - jl);
    % Loads DOFs
    loaddof = [3 * loadnid(:) - 2, 3 * loadnid(:) - 1, 3 * loadnid(:)];

    %% Characteristic matrix
    f_components.w = w;
    % Eigen decomposition
    [f_components.Q_original, f_components.SIGMA_original] = eig(f_components.w);

    f_components.T_original = f_components.Q_original * sqrt(f_components.SIGMA_original) * f_components.Q_original.';

    f_components.T_original_inverse = inv(f_components.T_original);

    f_components.T = f_components.Q_original.' * diag(1 ./ diag(sqrt(f_components.SIGMA_original)));

    % Num. of uncertain loads
    f_components.n = 1;
    %
    n_rand_f = 3 * f_components.n;
    % DOFs
    f_components.DOF = loaddof;
    % Nominal value
    f_components.f0 = [2.5; 0; 0];
    % Total DOFs
    f_components.total_DOFs = ndof;

    %% USER-DEFINED SUPPORT FIXED DOFs
    %% Cantilever beam fixed at left facet
    % [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
    % fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
    % fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs

    %% 3D wheel
    iif = [nelx nelx nelx nelx]; jf = [0 0 nely nely]; kf = [0 nelz 0 nelz];
    fixednid = kf * (nelx + 1) * (nely + 1) + iif * (nely + 1) + (nely + 1 - jf);
    fixeddof = [3 * fixednid(:); 3 * fixednid(:) - 1; 3 * fixednid(:) - 2];

    %%
    freedofs = setdiff(1:ndof, fixeddof);
    KE = lk_H8(nu);
    nodegrd = reshape(1:(nely + 1) * (nelx + 1), nely + 1, nelx + 1);
    nodeids = reshape(nodegrd(1:end - 1, 1:end - 1), nely * nelx, 1);
    nodeidz = 0:(nely + 1) * (nelx + 1):(nelz - 1) * (nely + 1) * (nelx + 1);
    nodeids = repmat(nodeids, size(nodeidz)) + repmat(nodeidz, size(nodeids));
    edofVec = 3 * nodeids(:) + 1;
    edofMat = repmat(edofVec, 1, 24) + ...
        repmat([0 1 2 3 * nely + [3 4 5 0 1 2] -3 -2 -1 ...
                            3 * (nely + 1) * (nelx + 1) + [0 1 2 3 * nely + [3 4 5 0 1 2] -3 -2 -1]], nele, 1);
    iK = reshape(kron(edofMat, ones(24, 1))', 24 * 24 * nele, 1);
    jK = reshape(kron(edofMat, ones(1, 24))', 24 * 24 * nele, 1);

    %% PREPARE FILTER
    iH = ones(nele * (2 * (ceil(rmin) - 1) + 1)^2, 1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;

    for k1 = 1:nelz

        for i1 = 1:nelx

            for j1 = 1:nely
                e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1;

                for k2 = max(k1 - (ceil(rmin) - 1), 1):min(k1 + (ceil(rmin) - 1), nelz)

                    for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)

                        for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                            e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2;
                            k = k + 1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0, rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2 + (k1 - k2)^2));
                        end

                    end

                end

            end

        end

    end

    H = sparse(iH, jH, sH);
    Hs = sum(H, 2);

    %% INITIALIZE ITERATION
    rho = repmat(volfrac, [nely, nelx, nelz]);
    rho_Phys = rho;
    iter = 0;
    change = 1;

    %% MMA parameters setup
    n_constraint = 1;
    n_design_variable = nelx * nely * nelz;
    rho_min = zeros(nely, nelx, nelz);
    rho_max = ones(nely, nelx, nelz);
    rho_old1_column = zeros(nelx * nely * nelz, 1);
    rho_old2_column = zeros(nelx * nely * nelz, 1);
    lower_asymptotes = zeros(nelx * nely * nelz, 1);
    upper_asymptotes = zeros(nelx * nely * nelz, 1);
    a0_MMA_input = 1;
    a_MMA_input = zeros(n_constraint, 1);
    c_MMA_input = 1000 * ones(n_constraint, 1);
    d_MMA_input = ones(n_constraint, 1);

    %%
    n_rand_var = n_rand_f + eole_order_e;

    %%
    mesh_info3D = element_xyz(nelx, nely, nelz);
    mesh_info3D.ndof = ndof;
    mesh_info3D.nele = nele;
    mesh_info3D.iK = iK;
    mesh_info3D.jK = jK;
    mesh_info3D.freedofs = freedofs;

    %% Material properties
    material_properties.edofMat = edofMat;
    material_properties.KE = KE;
    material_properties.Emin = Emin;
    material_properties.E0 = E0 * ones(nely, nelx, nelz);

    %% Generating Sparse Grid quadrature nodes and the corresponding weights
    [SG_nodes, SG_weights] = nwspgr('KPN', n_rand_var, k_level);

    %% Discretization of Young's modulus random field using EOLE
    if eole_order_e > 0
        [~, phi3D] = eole3D(mesh_info3D, eole_order_e, correlation_length_factor);
    end

    %%
    if switch_PCE == 1

        %% Generating PCE multi-index
        [PCE_dimension, PCE_multi_index] = generate_pce_multi_index(n_rand_var, pce_order);

        %% Calculating PCE normalization factors for compliance
        PCE_gamma = zeros(PCE_dimension, 1);

        for ii = 1:PCE_dimension
            PCE_gamma(ii, 1) = pce_normal_factor(PCE_multi_index(ii, :));
        end

        %% Calculating PCE base functions values at quadrature nodes
        PCE_base_nodal_value_SG = PCE_base_nodal(SG_nodes, PCE_dimension, PCE_multi_index);

        %% Initialization of PCE factors
        c_hat_SG = zeros(PCE_dimension, 1);
        Dc_hat_SG = zeros(nely, nelx, nelz, PCE_dimension);

    end

    %%
    if switch_UDR == 1
        %% Generating Gauss-Hermite quadrature nodes and the corresponding weights
        [GH_nodes, GH_weights] = gauss_hermite_node(n_tp_grid);

    end

    %% Iteration histories
    iter_history = zeros(max_iter, 6);
    rho_history = zeros(nely, nelx, nelz, max_iter + 1);
    rho_history(:, :, :, 1) = rho_Phys;

    %% START ITERATION
    while change > min_change && iter <= max_iter
        iter = iter + 1;
        % 	tic
        %% Increase penalty factor by 0.25 every 10 steps after the 90th step
        % if iter >= 30 && mod(iter, 5) == 0
        % penal = penal + 0.5;
        % end

        %% Creating function handle for FEA
        FEA_handle = @(xi_all) ...
            FEA3D(material_properties, mesh_info3D, rho_Phys, penal, f_components, n_rand_f, eole_order_e, err_magnitude_e, phi3D, xi_all);

        if switch_PCE == 1

            %% Calculating structural response at SPARSE GRID quadrature nodes
            [c_nodal_value_SG, Dc_nodal_value_SG] = f_integration_nodal3D(FEA_handle, SG_nodes, nelx, nely, nelz);

            %% Calculating PCE coefficients
            for ii = 1:PCE_dimension
                c_hat_SG(ii, 1) = (1 / PCE_gamma(ii)) * sum(c_nodal_value_SG .* PCE_base_nodal_value_SG(:, ii) .* SG_weights);

                Dc_nodal_value_SG_temp = Dc_nodal_value_SG;

                for jj = 1:size(SG_nodes, 1)
                    Dc_nodal_value_SG_temp(:, :, :, jj) = ...
                        Dc_nodal_value_SG_temp(:, :, :, jj) * PCE_base_nodal_value_SG(jj, ii) * SG_weights(jj, 1);
                end

                Dc_hat_SG(:, :, :, ii) = (1 / PCE_gamma(ii, 1)) * sum(Dc_nodal_value_SG_temp, 4);
            end

            %% Mean value and the mean value sensitivity of compliance (SPARSE GRID)
            c_mean_PCE_SG = c_hat_SG(1, 1);
            Dc_mean_PCE_SG = Dc_hat_SG(:, :, :, 1);

            %% STD and STD sensitivity of compliance (SPARSE GRID)
            c_STD_PCE_SG = sqrt(PCE_gamma(2:end, 1).' * (c_hat_SG(2:end).^2));

            Dc_STD_PCE_SG_temp = Dc_hat_SG;

            for ii = 1:PCE_dimension
                Dc_STD_PCE_SG_temp(:, :, :, ii) = PCE_gamma(ii, 1) * c_hat_SG(ii, 1) * Dc_hat_SG(:, :, :, ii);
            end

            Dc_STD_PCE_SG = (1 / c_STD_PCE_SG) * sum(Dc_STD_PCE_SG_temp, 4);

            %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            obj_PCE_SG = c_mean_PCE_SG + beta0 * c_STD_PCE_SG;
            Dobj_PCE_SG = Dc_mean_PCE_SG + beta0 * Dc_STD_PCE_SG;

            iter_history(iter, 1:2) = [c_mean_PCE_SG, c_STD_PCE_SG];

        end

        if switch_UDR == 1

            if n_MC_samples > 0
                [c_mean_UDR, Dc_mean_UDR, c_STD_UDR, Dc_STD_UDR, c_mean_MC, ~, c_STD_MC, ~] = ...
                    udr3D(mesh_info3D, n_rand_f, eole_order_e, GH_nodes, GH_weights, FEA_handle, n_MC_samples);

                iter_history(iter, 3:6) = [c_mean_UDR, c_STD_UDR, c_mean_MC, c_STD_MC];

            else

                [c_mean_UDR, Dc_mean_UDR, c_STD_UDR, Dc_STD_UDR] = ...
                    udr3D(mesh_info3D, n_rand_f, eole_order_e, GH_nodes, GH_weights, FEA_handle, n_MC_samples);

                iter_history(iter, 3:4) = [c_mean_UDR, c_STD_UDR];

            end

            obj_UDR = c_mean_UDR + beta0 * c_STD_UDR;
            Dobj_UDR = Dc_mean_UDR + beta0 * Dc_STD_UDR;

        end

        if switch_PCE == 1
            obj = obj_PCE_SG;
            Dobj = Dobj_PCE_SG;
        elseif switch_UDR == 1
            obj = obj_UDR;
            Dobj = Dobj_UDR;
        end

        %% ------------ Update design variables using the OC optimizer ------------ %%
        % %% FILTERING AND MODIFICATION OF SENSITIVITIES
        %         Dv = ones(nely, nelx, nelz);
        %         Dobj(:) = H * (Dobj(:) ./ Hs);
        %         Dv(:) = H * (Dv(:) ./ Hs);

        % %% OPTIMALITY CRITERIA UPDATE
        %         l1 = 0; l2 = 1e9; move = 0.2;

        %         while (l2 - l1) / (l1 + l2) > 1e-3
        %             lmid = 0.5 * (l2 + l1);
        %             rho_new = max(0, max(rho - move, min(1, min(rho + move, rho .* sqrt(-Dobj ./ Dv / lmid)))));
        %             rho_Phys(:) = (H * rho_new(:)) ./ Hs;

        %             if sum(rho_Phys(:)) > volfrac * nele
        %                 l1 = lmid;
        %             else
        %                 l2 = lmid;
        %             end

        %         end

        % change = max(abs(rho_new(:) - rho(:)) ./ rho(:));
        %         rho = rho_new;

        %%% ------------ End of OC updating block ------------ %%%

        %% ------------ Update design variables using MMA ------------ %%
        %% FILTERING AND MODIFICATION OF SENSITIVITIES
        Dobj(:) = H * (Dobj(:) ./ Hs);
        Dobj2 = 0 * Dobj;
        v = sum(sum(sum(rho_Phys))) - nelx * nely * nelz * volfrac;
        Dv = ones(nely, nelx, nelz);
        Dv(:) = H * (Dv(:) ./ Hs);
        Dv2 = 0 * Dv;

        

        rho_old = rho(:);

        [rho_new_column, ~, ~, ~, ~, ~, ~, ~, ~, lower_asymptotes, upper_asymptotes] = ...
            mmasub(...
            n_constraint, n_design_variable, iter, ...
            rho(:), rho_min(:), rho_max(:), rho_old1_column, rho_old2_column, ...
            obj, Dobj(:), Dobj2(:), ...
            v, Dv(:), Dv2(:), ...
            lower_asymptotes, upper_asymptotes, ...
            a0_MMA_input, a_MMA_input, c_MMA_input, d_MMA_input);

        rho = reshape(rho_new_column, nely, nelx, nelz);

        %% Apply density filter on the design variables vector
        rho_Phys = (H * rho_new_column) ./ Hs;
        rho_Phys = reshape(rho_Phys, nely, nelx, nelz);

        %% Update MMA parameters
        rho_old2_column = rho_old1_column;
        rho_old1_column = rho_old;

        change = max(abs(rho_new_column(:) - rho_old(:)) ./ rho(:));

        %%% ------------ End of MMA updating block ------------ %%%

        rho_history(:, :, :, iter + 1) = rho_Phys;

        %% PRINT RESULTS
        fprintf('\n It.:%5i  Vol.:%7.3f  ch.:%7.3f\n', iter, mean(rho_Phys(:)), change);
        fprintf('\nObjective = %.3f\n', obj);
        fprintf('\nMean_SG     STD_SG    Mean_UDR    STD_UDR    Mean_MC    STD_MC\n');
        disp(iter_history(iter, :));

        %% PLOT DENSITIES
        %     if displayflag == 1
        %         clf; display_3D(rho_Phys);
        %     end % #ok<UNRCH>

        %%
        save(cmd_output, 'rho_history', 'iter_history');

        %    toc
    end

    % fprintf('\nSolution converged after %d iteration steps.\n', iter);
    %
    % clf; display_3D(rho_Phys);

end % End of rto_3D()
