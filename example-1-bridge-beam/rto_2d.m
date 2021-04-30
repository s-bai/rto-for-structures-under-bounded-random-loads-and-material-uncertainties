function rto_2d(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid, n_MC_samples)
%% Robust topology optimization against uncertainties
% This code implements robust topology optimization of 2-dimensional
% structures against loads and material property uncertainties.

% Written by Song Bai.
% Supervised by Prof. Zhan Kang.

%   The "setprod.m" subroutine is developed by Mukhtar Ullah.
%   https://www.mathworks.com/matlabcentral/fileexchange/5898-setprod

%   The sparse grid and weight subroutine "nwspgr.m" is developed by Florian Heiss and Viktor Winschel.
%   http://www.sparse-grids.de/

% ALL VECTORS ARE STORED COLUMN-WISELY

% Input:
%   nelx:           Num. of elements in direction x;
%   nely:           Num. of elements in direction y;
%   volfrac:        Volume fraction restriction;
%   penalty:        Penalty factor;
%   rmin:           Filter radius;

%   k_level:        Accuracy level;
%   n_tp_grid:      Num. of 1-dimensional base nodes for tensor product grid
%   n_MC_samples:     Num. of Monte Carlo simulation samples for statistical moments

%   (Other input parameters are loaded from 'rto_2d.dat')

% Command example:
%   rto_2d(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid)
%   rto_2d( 100,   50,     0.5,       3,  1.5,       3,         5)


%% Switch for PCE (TENSOR PRODUCT and SPARSE GRID) and UDR methods (0 for no, 1 for yes)
switch_PCE = 1;
switch_TP = 1;
switch_SG = 1;
switch_UDR = 0;
switch_sensitivity_TP = 0;
switch_sensitivity_SG = 1;
switch_sensitivity_MC = 0;


%% Convergence criteria control parameters
max_iter_step = 200; % Maximum iteration steps
max_change = 0.01;


%% Start recording log file
clc;
cmd_file_name = 'rto_2d.log';
cmd_name = 'rto_2d';
cmd_parameters = [nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid];
diary (cmd_file_name);
cmd_output = cmd_history(cmd_name, cmd_parameters);

%% Load parameters from file
file_ID = fopen('rto_2d.dat', 'r'); % Open file
temp_file_read = textscan(file_ID, '%s %f', 'CommentStyle', '//');

temp = temp_file_read{2}(1:23);

%% Save parameters to a struct for simplicity
problem_parameters.nu = temp(1);

problem_parameters.e = temp(2);
e = problem_parameters.e*ones(nely + 1, nelx + 1);
problem_parameters.err_magnitude_e = temp(3);
problem_parameters.correlation_length_e = temp(4);

problem_parameters.t = temp(5);
t = problem_parameters.t*ones(nely + 1, nelx + 1);
problem_parameters.err_magnitude_t = temp(6);
problem_parameters.correlation_length_t = temp(7);

problem_parameters.eole_order_e = temp(8);
problem_parameters.eole_order_t = temp(9);
problem_parameters.pce_order = temp(10);
% problem_parameters.pce_order_t = temp(11);
problem_parameters.beta0 = temp(11);
problem_parameters.beta1 = temp(12);

% Characteristic matrix of the ellipsoid #1
temp1 = str2double(temp_file_read{1}(13:14));
temp2 = temp_file_read{2}(13:14);
ellips_char_parameter(1).w = horzcat(temp1, temp2);

% Parameters for the ellipsoid model #1
ellips_char_parameter(1).num_switch = 1;
ellips_char_parameter(1).fx_nominal = temp_file_read{2}(15);
ellips_char_parameter(1).fy_nominal = temp_file_read{2}(16);
ellips_char_parameter(1).size_value = -temp_file_read{2}(17);

% Temp variable for marking parameters of ellipoid model #2
load_switch = zeros(2, 1);
load_switch(1, 1) = temp_file_read{2}(18);

% Characteristic matrix of the ellipsoid #2
temp1 = str2double(temp_file_read{1}(19:20));
temp2 = temp_file_read{2}(19:20);
ellips_char_parameter(2).w = horzcat(temp1, temp2);

% Parameters for the ellipsoid model #2
ellips_char_parameter(2).num_switch = load_switch(1, 1);
ellips_char_parameter(2).fx_nominal = temp_file_read{2}(21);
ellips_char_parameter(2).fy_nominal = temp_file_read{2}(22);
ellips_char_parameter(2).size_value = -temp_file_read{2}(23);

% Temp variable for marking parameters of ellipoid model #3
load_switch(2, 1) = temp_file_read{2}(24);

% Characteristic matrix of the ellipsoid #3
temp1 = str2double(temp_file_read{1}(25:26));
temp2 = temp_file_read{2}(25:26);
ellips_char_parameter(3).w = horzcat(temp1, temp2);

% Parameters for the ellipsoid model #3
ellips_char_parameter(3).num_switch = load_switch(2, 1);
ellips_char_parameter(3).fx_nominal = temp_file_read{2}(27);
ellips_char_parameter(3).fy_nominal = temp_file_read{2}(28);
ellips_char_parameter(3).size_value = -temp_file_read{2}(29);

% Boundaries of uniform distribution
uniform_dist_boundary_e = zeros(1, 2);
uniform_dist_boundary_e(1) = temp_file_read{2}(30);
uniform_dist_boundary_e(2) = temp_file_read{2}(31);

uniform_dist_boundary_t = zeros(1, 2);
uniform_dist_boundary_t(1) = temp_file_read{2}(32);
uniform_dist_boundary_t(2) = temp_file_read{2}(33);


%%
fprintf('------------------------ Parameters setup ------------------------\n');
fprintf('\nNu = %.1f\n\nNominal Youngs modulus = %.1f, Error magnitude = %.1f\n', ...
    problem_parameters.nu, ...
    problem_parameters.e, ...
    problem_parameters.err_magnitude_e);

if problem_parameters.eole_order_e + problem_parameters.eole_order_t > 0
    fprintf('------------------------ Random field discretization ------------------------\n\n');
end

if problem_parameters.eole_order_t > 0
    fprintf('\nNominal thickness = %.1f, Error magnitude = %.1f\n', ...
        problem_parameters.t, ...
        problem_parameters.err_magnitude_t);
end

fprintf('\n\nYoungs modulus Correlation length = %.1f, Youngs modulus EOLE order = %d\n\n', ...
    problem_parameters.correlation_length_e, ...
    problem_parameters.eole_order_e);
fprintf('PCE order = %d\n\n', problem_parameters.pce_order);

if problem_parameters.eole_order_t > 0
    fprintf('\nThickness correlation length = %.1f, Thickness EOLE order = %d\n\n', ...
        problem_parameters.correlation_length_t, ...
        problem_parameters.eole_order_t);
end


%% ----------------- Boundary conditions -----------------
% All DOFs
allDOFs = 1:2*(nely + 1)*(nelx + 1);

%% -----------------  Cantilever beam -----------------
% boundary_condition.fixedDOFs = 1:1:2*(nely + 1);

%% ----------------- Simply supported beam (For MBB beam, set nelx/nely to 6:1) -----------------
boundary_condition.fixedDOFs = [1:2, 2*(nely + 1)*nelx + 2];

%% ----------------- Double-side hinged beam -----------------
% boundary_condition.fixedDOFs = [2*(nely + 1) - 1, 2*(nely + 1), 2*(nelx + 1)*(nely + 1) - 1, 2*(nelx + 1)*(nely + 1)];

%% ----------------- Two side fixed beam -----------------
% boundary_condition.fixedDOFs = [1:1:2*(nely + 1), 2*(nely + 1)*nelx + 1:2*(nelx + 1)*(nely + 1)];

%% ----------------- Double side mid hinged beam -----------------
% boundary_condition.fixedDOFs = [nely + 1, nely + 2, 2*(nely + 1)*nelx + nely + 1, 2*(nely + 1)*nelx + nely + 2];

%% ----------------- Rectangular plate with lower side fixed -----------------
% temp_node_num = linspace(nely + 1, (nelx + 1)*(nely + 1), nelx + 1);
% boundary_condition.fixedDOFs = sort([2*temp_node_num - 1, 2*temp_node_num]);

%% ----------------- Free DOFs -----------------
boundary_condition.freeDOFs = setdiff(allDOFs, boundary_condition.fixedDOFs);

%% ----------------- Nominal loads -----------------
%   The first column stores the DOFs of the loads
%   The second column stores the values
f_nominal = zeros(2, 2, sum(load_switch) + 1);

%% ----------------- Loads for cantilever beam with one loads -----------------
% % Loads at right mid-point
% f_nominal(:, :, 1) = [
%     2*nelx*(nely + 1) + nely + 1,  ellips_char_parameter(1).fx_nominal;
%     2*nelx*(nely + 1) + nely + 2,  ellips_char_parameter(1).fy_nominal];

%% ----------------- Loads for cantilever beam with two loads at the bottom -----------------
% % Loads at mid bottom point
% f_nominal(:, :, 1) = [
%     nelx*(nely + 1) + 1,  ellips_char_parameter(1).fx_nominal;
%     nelx*(nely + 1) + 2,  ellips_char_parameter(1).fy_nominal];
%
% % Loads at right bottom corner
% f_nominal(:, :, 2) = [
%     2*nelx*(nely + 1) + 1,  ellips_char_parameter(2).fx_nominal;
%     2*nelx*(nely + 1) + 2,  ellips_char_parameter(2).fy_nominal];

%% ----------------- Loads for MBB beam -----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     (nely + 1)*nelx/2 + 1,  ellips_char_parameter(1).fx_nominal;
%     (nely + 1)*nelx/2 + 2,  ellips_char_parameter(1).fy_nominal];
% %
% % % Loads #2
% f_nominal(:, :, 2) = [
%     (nely + 1)*nelx + 1,  ellips_char_parameter(2).fx_nominal;
%     (nely + 1)*nelx + 2,  ellips_char_parameter(2).fy_nominal];
% %
% % % Loads #3
% f_nominal(:, :, 3) = [
%     (nely + 1)*3/2*nelx + 1,  ellips_char_parameter(3).fx_nominal;
%     (nely + 1)*3/2*nelx + 2,  ellips_char_parameter(3).fy_nominal];

%% ----------------- Loads for double-side hinged beam -----------------
% Loads #1
f_nominal(:, :, 1) = [
    2*(nely + 1)*(nelx/2 + 1) - 1,  ellips_char_parameter(1).fx_nominal;
    2*(nely + 1)*(nelx/2 + 1),  ellips_char_parameter(1).fy_nominal];

%% ----------------- Loads for double-side fixed beam (Loads at 1/4 1/2 3/4)-----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     2*(nely + 1)*(nelx/4 + 1) - 1,  ellips_char_parameter(1).fx_nominal;
%     2*(nely + 1)*(nelx/4 + 1),  ellips_char_parameter(1).fy_nominal];
% %
% % % Loads #2
% f_nominal(:, :, 2) = [
%     2*(nely + 1)*(nelx/2 + 1) - 1,  ellips_char_parameter(2).fx_nominal;
%     2*(nely + 1)*(nelx/2 + 1),  ellips_char_parameter(2).fy_nominal];
% %
% % % Loads #3
% f_nominal(:, :, 3) = [
%     2*(nely + 1)*(nelx*3/4 + 1) - 1,  ellips_char_parameter(3).fx_nominal;
%     2*(nely + 1)*(nelx*3/4 + 1),  ellips_char_parameter(3).fy_nominal];

%% ----------------- Loads for double-side fixed beam (Loads at 1/3 1/2 2/3)-----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     2*(nely + 1)*(nelx/3 + 1) - 1,  ellips_char_parameter(1).fx_nominal;
%     2*(nely + 1)*(nelx/3 + 1),  ellips_char_parameter(1).fy_nominal];
% %
% % % Loads #2
% f_nominal(:, :, 2) = [
%     2*(nely + 1)*(nelx/2 + 1) - 1,  ellips_char_parameter(2).fx_nominal;
%     2*(nely + 1)*(nelx/2 + 1),  ellips_char_parameter(2).fy_nominal];
% %
% % % Loads #3
% f_nominal(:, :, 3) = [
%     2*(nely + 1)*(nelx*2/3 + 1) - 1,  ellips_char_parameter(3).fx_nominal;
%     2*(nely + 1)*(nelx*2/3 + 1),  ellips_char_parameter(3).fy_nominal];

%% ----------------- Loads for double-side fixed beam (Loads at 1/4 3/4)-----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     (nely + 1)*nelx/2 + nely + 1,  ellips_char_parameter(1).fx_nominal;
%     (nely + 1)*nelx/2 + nely + 2,  ellips_char_parameter(1).fy_nominal];
% %
% % % Loads #2
% f_nominal(:, :, 2) = [
%     (nely + 1)*3*nelx/2 + nely + 1,  ellips_char_parameter(2).fx_nominal;
%     (nely + 1)*3*nelx/2 + nely + 2,  ellips_char_parameter(2).fy_nominal];

%% ----------------- Loads for double-side hinge supported beam with loads at center -----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     nelx*(nely + 1) + nely + 1,  ellips_char_parameter(1).fx_nominal;
%     nelx*(nely + 1) + nely + 2,  ellips_char_parameter(1).fy_nominal];

%% ----------------- Loads for cantilever beam with loads at right upper and lower corners -----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     2*(nely + 1)*nelx + 1,  ellips_char_parameter(1).fx_nominal;
%     2*(nely + 1)*nelx + 2,  ellips_char_parameter(1).fy_nominal];
% %
% % % Loads #2
% f_nominal(:, :, 2) = [
%     2*(nely + 1)*(nelx + 1) - 1,  ellips_char_parameter(2).fx_nominal;
%     2*(nely + 1)*(nelx + 1),  ellips_char_parameter(2).fy_nominal];

%% ----------------- Loads for rectangular plate with lower side fixed -----------------
% % Loads #1
% f_nominal(:, :, 1) = [
%     (nely + 1)*nelx + 1,  ellips_char_parameter(1).fx_nominal;
%     (nely + 1)*nelx + 2,  ellips_char_parameter(1).fy_nominal];


%%
fclose(file_ID); % Close file

%% Num. of ellipsoid models
n_ellipse = 0;
for ii = 1:3
    n_ellipse = n_ellipse + ellips_char_parameter(ii).num_switch;
end
problem_parameters.n_ellipse = n_ellipse;

%% Num. of uncertain loads (Can be overriden to be zero)
rand_numbr_f = 2*n_ellipse;

%% Building loads component struct
f_components.n = rand_numbr_f;
f_components.DOF = zeros(2*n_ellipse, 1);
f_components.nominal_value = zeros(2*n_ellipse, 1);

for ii = 1:n_ellipse
    f_components.DOF(2*ii - 1:2*ii, 1) = f_nominal(:, 1, ii);
    f_components.nominal_value(2*ii - 1:2*ii, 1) = f_nominal(:, 2, ii);
end



%% Meshing
mesh_info = Meshing(nelx, nely);


%% Eigenvalue and eigenvector decomposition of characteristic matrix w
w = struct([]);
for ii = 1:n_ellipse
    w(ii).w_matrix = ellips_char_parameter(ii).w;
    [w(ii).Q, w(ii).SIGMA] = eig(w(ii).w_matrix);
    w(ii).T = w(ii).Q*sqrt(w(ii).SIGMA)*w(ii).Q.';
    w(ii).T_inverse = inv(w(ii).T);
end


%% Calculating total EOLE orders and PCE orders
n_rand_var = rand_numbr_f + problem_parameters.eole_order_e + problem_parameters.eole_order_t;
n_pce_order = problem_parameters.pce_order;


%% Generating Gauss-Hermite quadrature nodes and the corresponding weights
[GH_nodes, GH_weights] = gauss_hermite_node(n_tp_grid);

%%
if switch_PCE == 1
    %% Generating Sparse grid nodes and tensor product nodes for tensor product and sparse grid methods
    [TP_nodes, TP_weights] = tp_nodes_generation(n_rand_var, n_tp_grid, GH_nodes, GH_weights);
    
    %% Generating Sparse Grid quadrature nodes and the corresponding weights
    [SG_nodes, SG_weights] = nwspgr('KPN', n_rand_var, k_level);
    
    %%
    if problem_parameters.eole_order_t > 0
        %% Generating SG nodes and TP nodes for volume calculation
        [TP_nodes_vol, TP_weights_vol] = tp_nodes_generation(problem_parameters.eole_order_t, n_tp_grid, GH_nodes, GH_weights);
        
        %% Generating Sparse Grid quadrature nodes and the corresponding weights
        [SG_nodes_vol, SG_weights_vol] = nwspgr('KPN', problem_parameters.eole_order_t, k_level);
        
    end
    
end



%% Calculating mode functions of the Young's modulus random field
if problem_parameters.eole_order_e > 0
    [phi_e] = eole_mode_fun(mesh_info, problem_parameters.correlation_length_e, problem_parameters.eole_order_e);
else
    phi_e = 0;
end


%% Calculating mode functions of the thickness random field
if problem_parameters.eole_order_t > 0
    [phi_t] = eole_mode_fun(mesh_info, problem_parameters.correlation_length_t, problem_parameters.eole_order_t);
else
    phi_t = 0;
end


%% Generating PCE multi-index
if switch_PCE == 1
    [compliance_pce_dimension, compliance_pce_multi_index] = generate_pce_multi_index(n_rand_var, n_pce_order);
    
    %% Calculating PCE normalization factors for compliance
    pce_gamma_compliance = zeros(compliance_pce_dimension, 1);
    
    for ii = 1:compliance_pce_dimension
        pce_gamma_compliance(ii, 1) = pce_normal_factor(compliance_pce_multi_index(ii, :));
    end
    
    if problem_parameters.eole_order_t > 0
        [volume_pce_dimension, volume_pce_multi_index] = ...
            generate_pce_multi_index(problem_parameters.eole_order_t, problem_parameters.pce_order);
        
        %% Calculating PCE normalization factors for volume
        pce_gamma_volume = zeros(volume_pce_dimension, 1);
        
        for ii = 1:volume_pce_dimension
            pce_gamma_volume(ii, 1) = pce_normal_factor(volume_pce_multi_index(ii, :));
        end
    end
end



%% Prepare filter (From 'top88' by E. Andreassen et al.)
iH = ones(nelx*nely*(2*(ceil(rmin) - 1) + 1)^2, 1);
jH = ones(size(iH));
sH = zeros(size(iH));
k_filter = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely + j1;
        for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)
            for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                e2 = (i2 - 1)*nely + j2;
                k_filter = k_filter + 1;
                iH(k_filter) = e1;
                jH(k_filter) = e2;
                sH(k_filter) = max(0, rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2));
            end
        end
    end
end
H = sparse(iH, jH, sH);
Hs = sum(H, 2);


%%
if switch_PCE == 1
    %% Initialization of PCE factors
    if switch_TP == 1
        c_hat_TP = zeros(compliance_pce_dimension, 1);
        Dc_hat_TP = zeros(nely, nelx, compliance_pce_dimension);
    end
    
    if switch_SG == 1
        c_hat_SG = zeros(compliance_pce_dimension, 1);
        Dc_hat_SG = zeros(nely, nelx, compliance_pce_dimension);
    end
    
    %%
    if problem_parameters.eole_order_t > 0
        if switch_TP == 1
            v_hat_TP = zeros(volume_pce_dimension, 1);
            Dv_hat_TP = zeros(nely, nelx, volume_pce_dimension);
        end
        
        if switch_SG == 1
            v_hat_SG = zeros(volume_pce_dimension, 1);
            Dv_hat_SG = zeros(nely, nelx, volume_pce_dimension);
        end
        
    end
    
end



%% Initialization of optimization parameters
rho_initial = volfrac;
rho = repmat(rho_initial, nely, nelx); % Element density matrix

iter = 1; % Iteration step
change = 1; % Convergency criteria

%% [c_mean_PCE_TP, c_STD_PCE_TP, c_mean_PCE_SG, c_STD_PCE_SG, c_mean_UDR, c_STD_UDR, c_mean_MC, c_STD_MC]
objective_iter_history = zeros(max_iter_step, 8);

%% [v_mean_frac_PCE_TP, v_STD_frac_PCE_TP, v_mean_frac_PCE_SG, v_STD_frac_PCE_SG, v_mean_fraction_UDR, v_STD_fraction_UDR]
constraint_iter_history = zeros(max_iter_step, 6);

%%
rho_history = zeros(nely, nelx, max_iter_step + 1);
rho_history(:, :, 1) = rho; % Initial design as the first step

%% Parameters for MMA
n_constraint = 1;
n_design_variable = nelx*nely;
rho_min = zeros(nely, nelx);
rho_max = ones(nely, nelx);
rho_old1_column = zeros(nelx*nely, 1);
rho_old2_column = zeros(nelx*nely, 1);
lower_asymptotes = zeros(nelx*nely, 1);
upper_asymptotes = zeros(nelx*nely, 1);
a0_MMA_input = 1;
a_MMA_input = zeros(n_constraint, 1);
c_MMA_input = 1000*ones(n_constraint, 1);
d_MMA_input = ones(n_constraint, 1);



%% Creating function handle array for PCE base functions
if switch_PCE == 1
    if switch_TP == 1
        FEA_PCE_base_nodal_value_TP = PCE_base_nodal(TP_nodes, compliance_pce_dimension, compliance_pce_multi_index);
        
        if problem_parameters.eole_order_t > 0
            Vol_PCE_base_nodal_value_TP = PCE_base_nodal(TP_nodes_vol, volume_pce_dimension, volume_pce_multi_index);
        end
    end
    
    if switch_SG == 1
        FEA_PCE_base_nodal_value_SG = PCE_base_nodal(SG_nodes, compliance_pce_dimension, compliance_pce_multi_index);
        
        if problem_parameters.eole_order_t > 0
            Vol_PCE_base_nodal_value_SG = PCE_base_nodal(SG_nodes_vol, volume_pce_dimension, volume_pce_multi_index);
        end
    end
    
end





%% Start optimization
figure('Name', 'Topology layout evolution', 'NumberTitle', 'off');

while (change > max_change) && (iter <= max_iter_step)
    
    fprintf('\n------------------------ Iter. %d ------------------------\n', iter);
    
    %% Increase penalty factor by 0.2 every 10 steps after the 100th step
    if iter > 100 && mod(iter, 10) == 0
        penalty = penalty + 0.25;
    end
    
    
    %% Creating function handle for FEA
    FEA_handle = @(xi_all) FEA(mesh_info, e, problem_parameters.eole_order_e, ...
        t, problem_parameters.eole_order_t, ...
        problem_parameters.nu, boundary_condition, ...
        f_components, rand_numbr_f, rho, penalty, w, ...
        phi_e, problem_parameters.err_magnitude_e, uniform_dist_boundary_e, ...
        phi_t, problem_parameters.err_magnitude_t, uniform_dist_boundary_t, xi_all);
    
    %%
    if problem_parameters.eole_order_t > 0
        %% Creating functon handle for volume
        volume_handle = @(xi_all) total_volume(mesh_info, t, rho, ...
            problem_parameters.eole_order_t, xi_all, phi_t, problem_parameters.err_magnitude_t, uniform_dist_boundary_t);
    end
    
  
    
    %% Calculating structural responses at quadrature nodes by PCE
    if switch_PCE == 1
        if switch_TP == 1
            %% Calculating structural response at TENSOR PRODUCT quadrature nodes
            [c_nodal_value_TP, Dc_nodal_value_TP] = f_integration_nodal(FEA_handle, TP_nodes, nelx, nely);
            
            for ii = 1:compliance_pce_dimension
                c_hat_TP(ii, 1) = (1/pce_gamma_compliance(ii, 1))*sum(c_nodal_value_TP.*FEA_PCE_base_nodal_value_TP(:, ii).*TP_weights);
                
                Dc_nodal_value_TP_temp = Dc_nodal_value_TP;
                for jj = 1:size(TP_nodes, 1)
                    Dc_nodal_value_TP_temp(:, :, jj) = ...
                        Dc_nodal_value_TP_temp(:, :, jj)*FEA_PCE_base_nodal_value_TP(jj, ii)*TP_weights(jj, 1);
                end
                Dc_hat_TP(:, :, ii) = (1/pce_gamma_compliance(ii, 1))*sum(Dc_nodal_value_TP_temp, 3);
                
            end
            
            %% Mean value and mean value sensitivity of compliance (TENSOR PRODUCT)
            c_mean_PCE_TP = c_hat_TP(1, 1);
            Dc_mean_PCE_TP = Dc_hat_TP(:, :, 1);
            
            %% STD and STD sensitivity of compliance (TENSOR PRODUCT)
            c_STD_PCE_TP = sqrt(pce_gamma_compliance(2:end, 1).'*(c_hat_TP(2:end).^2));
            
            Dc_STD_PCE_TP_temp = Dc_hat_TP;
            for ii = 1:compliance_pce_dimension
                Dc_STD_PCE_TP_temp(:, :, ii) = pce_gamma_compliance(ii, 1)*c_hat_TP(ii, 1)*Dc_hat_TP(:, :, ii);
            end
            Dc_STD_PCE_TP = (1/c_STD_PCE_TP)*sum(Dc_STD_PCE_TP_temp, 3);
            
            objective_iter_history(iter, 1:2) = [c_mean_PCE_TP, c_STD_PCE_TP];
        end
        
        if switch_SG == 1
            %% Calculating structural response at SPARSE GRID quadrature nodes
            [c_nodal_value_SG, Dc_nodal_value_SG] = f_integration_nodal(FEA_handle, SG_nodes, nelx, nely);
            
%             c_nodal_value_SG = 100*c_nodal_value_SG;
%             Dc_nodal_value_SG = 100*Dc_nodal_value_SG;
            
            for ii = 1:compliance_pce_dimension
                c_hat_SG(ii, 1) = (1/pce_gamma_compliance(ii))*sum(c_nodal_value_SG.*FEA_PCE_base_nodal_value_SG(:, ii).*SG_weights);
                
                Dc_nodal_value_SG_temp = Dc_nodal_value_SG;
                for jj = 1:size(SG_nodes, 1)
                    Dc_nodal_value_SG_temp(:, :, jj) = ...
                        Dc_nodal_value_SG_temp(:, :, jj)*FEA_PCE_base_nodal_value_SG(jj, ii)*SG_weights(jj, 1);
                end
                Dc_hat_SG(:, :, ii) = (1/pce_gamma_compliance(ii, 1))*sum(Dc_nodal_value_SG_temp, 3);
                
            end
            
            %% Mean value and the mean value sensitivity of compliance (SPARSE GRID)
            c_mean_PCE_SG = c_hat_SG(1, 1);
            Dc_mean_PCE_SG = Dc_hat_SG(:, :, 1);
            
            %% STD and STD sensitivity of compliance (SPARSE GRID)
            c_STD_PCE_SG = sqrt(pce_gamma_compliance(2:end, 1).'*(c_hat_SG(2:end).^2));
            
            Dc_STD_PCE_SG_temp = Dc_hat_SG;
            for ii = 1:compliance_pce_dimension
                Dc_STD_PCE_SG_temp(:, :, ii) = pce_gamma_compliance(ii, 1)*c_hat_SG(ii, 1)*Dc_hat_SG(:, :, ii);
            end
            Dc_STD_PCE_SG = (1/c_STD_PCE_SG)*sum(Dc_STD_PCE_SG_temp, 3);
            
            objective_iter_history(iter, 3:4) = [c_mean_PCE_SG, c_STD_PCE_SG];
        end
        
    end
    
    
    %% Calculating structural responses at quadrature nodes by UDR w/o Monte Carlo simulation
    if switch_UDR == 1
        if n_MC_samples > 0
            %% Implementing both UDR and Monte Carlo simulation
            [c_mean_UDR, Dc_mean_UDR, c_STD_UDR, Dc_STD_UDR, c_mean_MC, Dc_mean_MC, c_STD_MC, Dc_STD_MC] = ...
                udr(mesh_info, rand_numbr_f, problem_parameters.eole_order_e, problem_parameters.eole_order_t, ...
                GH_nodes, GH_weights, FEA_handle, n_MC_samples);
                        
            objective_iter_history(iter, 5:8) = [c_mean_UDR, c_STD_UDR, c_mean_MC, c_STD_MC];
            
        else
            
            %% Implementing UDR only
            [c_mean_UDR, Dc_mean_UDR, c_STD_UDR, Dc_STD_UDR] = ...
                udr(mesh_info, rand_numbr_f, problem_parameters.eole_order_e, problem_parameters.eole_order_t, ...
                GH_nodes, GH_weights, FEA_handle, n_MC_samples);
 
            objective_iter_history(iter, 5:6) = [c_mean_UDR, c_STD_UDR];
            
        end
    end
    
    
    %% Objective
    if switch_sensitivity_MC == 1
        %% Use the objective value and its sensitivity obtained by Monte Carlo
        rto_obj = c_mean_MC + problem_parameters.beta0*c_STD_MC;
        Drto_obj = Dc_mean_MC + problem_parameters.beta0*Dc_STD_MC;
        Drto_obj2 = 0*Drto_obj;
        
    elseif switch_sensitivity_TP == 1
        %% Use the objective value and its sensitivity obtained by PCE TENSOR PRODUCT
        rto_obj = c_mean_PCE_TP + problem_parameters.beta0*c_STD_PCE_TP;
        Drto_obj = Dc_mean_PCE_TP + problem_parameters.beta0*Dc_STD_PCE_TP;
        Drto_obj2 = 0*Drto_obj;
        
    elseif switch_sensitivity_SG == 1
        %% Use the sensitivity obtained by PCE SPARSE GRID
        rto_obj = c_mean_PCE_SG + problem_parameters.beta0*c_STD_PCE_SG;
        Drto_obj = Dc_mean_PCE_SG + problem_parameters.beta0*Dc_STD_PCE_SG;
        Drto_obj2 = 0*Drto_obj;
        
    else
        %% Use the sensitivity obtained by UDR
        rto_obj = c_mean_UDR + problem_parameters.beta0*c_STD_UDR;
        Drto_obj = Dc_mean_UDR + problem_parameters.beta0*Dc_STD_UDR;
        Drto_obj2 = 0*Drto_obj;
        
    end
    
    
    %% Calculating volume response
    if problem_parameters.eole_order_t == 0
        %% Ignoring thickness uncertainties
        rto_constraint = sum(sum(rho)) - nelx*nely*volfrac;
        Drto_constraint = ones(nely, nelx);
        Drto_constraint2 = 0*Drto_constraint;
        
        v_mean_fraction = sum(sum(rho))/(nelx*nely*volfrac);
        
    else
        %% Considering thickness uncertainties
        if switch_PCE == 1
            %% Calculating volume mean and std by PCE
            if switch_TP == 1
                %% Calculating volume at TENSOR PRODUCT nodes
                [v_nodal_value_TP, Dv_nodal_value_TP] = f_integration_nodal(volume_handle, TP_nodes_vol, nelx, nely);
                
                for ii = 1:volume_pce_dimension
                    v_hat_TP(ii, 1) = (1/pce_gamma_volume(ii, 1))*sum(v_nodal_value_TP.*Vol_PCE_base_nodal_value_TP(:, ii).*TP_weights_vol);
                    
                    Dv_nodal_value_TP_temp = Dv_nodal_value_TP;
                    for jj = 1:size(TP_nodes_vol, 1)
                        Dv_nodal_value_TP_temp(:, :, jj) = ...
                            Dv_nodal_value_TP_temp(:, :, jj)*Vol_PCE_base_nodal_value_TP(jj, ii)*TP_weights_vol(jj, 1);
                    end
                    Dv_hat_TP(:, :, ii) = (1/pce_gamma_volume(ii, 1))*sum(Dv_nodal_value_TP_temp, 3);
                end
                
                %% Mean value and mean value sensitivity of volume (TENSOR PRODUCT)
                v_mean_PCE_TP = v_hat_TP(1, 1);
                Dv_mean_PCE_TP = Dv_hat_TP(:, :, 1);
                
                %% STD and STD sensitivity of volume (TENSOR PRODUCT)
                v_STD_PCE_TP = sqrt(pce_gamma_volume(2:end, 1).'*(v_hat_TP(2:end).^2));
                
                Dv_STD_PCE_TP_temp = Dv_hat_TP;
                for ii = 1:volume_pce_dimension
                    Dv_STD_PCE_TP_temp(:, :, ii) = pce_gamma_volume(ii, 1)*v_hat_TP(ii, 1)*Dv_hat_TP(:, :, ii);
                end
                Dv_STD_PCE_TP = (1/v_STD_PCE_TP)*sum(Dv_STD_PCE_TP_temp, 3);
                
            end
            
            if switch_SG == 1
                %% Calculating volume at SPARSE GRID nodes
                [v_nodal_value_SG, Dv_nodal_value_SG] = f_integration_nodal(volume_handle, SG_nodes_vol, nelx, nely);
                
                for ii = 1:volume_pce_dimension
                    v_hat_SG(ii, 1) = (1/pce_gamma_volume(ii, 1))*sum(v_nodal_value_SG.*Vol_PCE_base_nodal_value_SG(:, ii).*SG_weights_vol);
                    
                    Dv_nodal_value_SG_temp = Dv_nodal_value_SG;
                    for jj = 1:size(SG_nodes_vol, 1)
                        Dv_nodal_value_SG_temp(:, :, jj) = ...
                            Dv_nodal_value_SG_temp(:, :, jj)*Vol_PCE_base_nodal_value_SG(jj, ii)*SG_weights_vol(jj, 1);
                    end
                    Dv_hat_SG(:, :, ii) = (1/pce_gamma_volume(ii, 1))*sum(Dv_nodal_value_SG_temp, 3);
                end
                
                %% Mean value and mean value sensitivity of volume (SPARSE GRID)
                v_mean_PCE_SG = v_hat_SG(1, 1);
                Dv_mean_PCE_SG = Dv_hat_SG(:, :, 1);
                
                %% STD and STD sensitivity of volume (SPARSE GRID)
                v_STD_PCE_SG = sqrt(pce_gamma_volume(2:end, 1).'*(v_hat_SG(2:end).^2));
                
                Dv_STD_PCE_SG_temp = Dv_hat_SG;
                for ii = 1:volume_pce_dimension
                    Dv_STD_PCE_SG_temp(:, :, ii) = pce_gamma_volume(ii, 1)*v_hat_SG(ii, 1)*Dv_hat_SG(:, :, ii);
                end
                Dv_STD_PCE_SG = (1/v_STD_PCE_SG)*sum(Dv_STD_PCE_SG_temp, 3);
                
            end
        end
        
        if switch_UDR == 1
            %% Calculating volume mean and std by UDR and Monte Carlo
            if n_MC_samples > 0
                [v_mean_UDR, Dv_mean_UDR, v_STD_UDR, Dv_STD_UDR, v_mean_MC, ~, v_STD_MC, ~] = ...
                    udr(mesh_info, 0, 0, ...
                    problem_parameters.eole_order_t, ...
                    GH_nodes, GH_weights, volume_handle, n_MC_samples);
                
                fprintf('\nVolume mean by UDR = %.2f, STD by UDR = %.2f\n', ...
                    v_mean_UDR/(volfrac*nelx*nely), v_STD_UDR/(volfrac*nelx*nely));
                
                fprintf('\nVolume mean by Monte Carlo = %.2f, STD by Monte Carlo = %.2f\n', ...
                    v_mean_MC/(volfrac*nelx*nely), v_STD_MC/(volfrac*nelx*nely));
                
            else
                [v_mean_UDR, Dv_mean_UDR, v_STD_UDR, Dv_STD_UDR] = ...
                    udr(mesh_info, 0, 0, ...
                    problem_parameters.eole_order_t, ...
                    GH_nodes, GH_weights, volume_handle, n_MC_samples);
                
                fprintf('\nVolume mean by UDR = %.2f, STD by UDR = %.2f\n', ...
                    v_mean_UDR/(volfrac*nelx*nely), v_STD_UDR/(volfrac*nelx*nely));
                
            end
            
        end
        
        %% Calculating constraints
        if switch_sensitivity_TP == 1
            %% Use the constraint value and its sensitivity obtained by PCE TENSOR PRODUCT
            rto_constraint = v_mean_PCE_TP + 1*v_STD_PCE_TP - volfrac*nelx*nely;
            Drto_constraint = Dv_mean_PCE_TP + 1*Dv_STD_PCE_TP;
            Drto_constraint2 = 0*Drto_constraint;
            
            v_mean_fraction_TP = (v_mean_PCE_TP)/(volfrac*nelx*nely);
            v_STD_fraction_TP = (1*v_STD_PCE_TP)/(volfrac*nelx*nely);
            
            constraint_iter_history(iter, 1:2) = [v_mean_fraction_TP, v_STD_fraction_TP];
            
        elseif switch_sensitivity_SG == 1
            %% Use the sensitivity obtained by PCE SPARSE GRID
            rto_constraint = v_mean_PCE_SG + 1*v_STD_PCE_SG - volfrac*nelx*nely;
            Drto_constraint = Dv_mean_PCE_SG + 1*Dv_STD_PCE_SG;
            Drto_constraint2 = 0*Drto_constraint;
            
            v_mean_fraction_SG = (v_mean_PCE_SG)/(volfrac*nelx*nely);
            v_STD_fraction_SG = (1*v_STD_PCE_SG)/(volfrac*nelx*nely);
            
            constraint_iter_history(iter, 3:4) = [v_mean_fraction_SG, v_STD_fraction_SG];
            
        else
            %% Use the sensitivity obtained by UDR
            rto_constraint = v_mean_UDR + 1*v_STD_UDR - volfrac*nelx*nely;
            Drto_constraint = Dv_mean_UDR + 1*Dv_STD_UDR;
            Drto_constraint2 = 0*Drto_constraint;
            
            v_mean_fraction_UDR = (v_mean_UDR)/(volfrac*nelx*nely);
            v_STD_fraction_UDR = (1*v_STD_UDR)/(volfrac*nelx*nely);
            
            constraint_iter_history(iter, 5:6) = [v_mean_fraction_UDR, v_STD_fraction_UDR];
            
        end
        
        
    end
    
    
    
    %% Objective sensitivity filtering
    Drto_obj(:) = H*(rho(:).*Drto_obj(:))./Hs./max(1e-3, rho(:));
    
    
    %% Updating design variables using MMA
    rho_old = rho(:);
    
    [rho_new_column, ~, ~, ~, ~, ~, ~, ~, ~, lower_asymptotes, upper_asymptotes] = ...
        mmasub(...
        n_constraint, n_design_variable, iter, ...
        rho(:), rho_min(:), rho_max(:), rho_old1_column, rho_old2_column, ...
        rto_obj, Drto_obj(:), Drto_obj2(:), ...
        rto_constraint, Drto_constraint(:), Drto_constraint2(:), ...
        lower_asymptotes, upper_asymptotes, ...
        a0_MMA_input, a_MMA_input, c_MMA_input, d_MMA_input);
    
    rho = reshape(rho_new_column, nely, nelx);
    rho_history(:, :, iter + 1) = rho;
    
    save(cmd_output, 'rho_history', 'objective_iter_history', 'constraint_iter_history');
    
    %% Print result
    if problem_parameters.eole_order_t == 0
        fprintf('\nc_mean_TP, c_STD_TP, c_mean_SG, c_STD_SG, c_mean_UDR, c_STD_UDR, c_mean_MC, c_STD_MC\n');
        disp(objective_iter_history(iter, :));
        
        fprintf('\nObj: %.3f Volume: %.3f change: %.3f\n', rto_obj, v_mean_fraction, change);
        
    else
        fprintf('\nc_mean_TP, c_STD_TP, c_mean_SG, c_STD_SG, c_mean_UDR, c_STD_UDR, c_mean_MC, c_STD_MC\n');
        disp(objective_iter_history(iter, :));
        
        fprintf('\nObjective: %.4f change: %.3f\n', rto_obj, change);
        
        fprintf('\nv_mean_TP, v_STD_TP, v_mean_SG, v_STD_SG, v_mean_UDR, v_STD_UDR]\n');
        disp(constraint_iter_history(iter, :));
        
    end
    
    %% Plotting topology layout (density)
    colormap(gray);
    imagesc(1 - rho);
    caxis([0 1]);
    axis equal;
    axis off;
    drawnow;
    
    rho_old2_column = rho_old1_column;
    rho_old1_column = rho_old;
    
    change = max(abs(rho(:) - rho_old(:)));
    
    iter = iter + 1;
end


%% Save iteration data to MAT file
iter = iter - 1;
fprintf('\nSolution converged after %d iteration steps.\n', iter);

save(cmd_output, 'problem_parameters', 'iter', 'ellips_char_parameter', 'n_ellipse', 'rand_numbr_f', ...
    'f_components', 'w', ...
    'boundary_condition', 'uniform_dist_boundary_e', 'uniform_dist_boundary_t', '-append');

%% End recording log file
diary off;

end % End of main function









