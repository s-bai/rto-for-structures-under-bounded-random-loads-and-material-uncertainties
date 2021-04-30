function [u, c, d_c] = FEA(mesh_info, e, eole_order_e, t, eole_order_t, nu, ...
        boundary_condition, f_components, rand_numbr_f, rho, penalty, ...
        w, ...
        phi_e, error_magnitude_e, uniform_dist_boundary_e, ...
        phi_t, error_magnitude_t, uniform_dist_boundary_t, xi_all)
    %% Finite element analysis subroutine

    % Input:
    %   mesh_info
    %       mesh_info.nelx: Num. of elements in direction x;
    %       mesh_info.nely: Num. of elements in direction y;
    %       mesh_info.ik: The row index vectors
    %       mesh_info.jk: The column index vectors
    %       mesh_info.element_node_num: Element node numbers;

    %   e: Young's modulus field;
    %   eole_order_e: EOLE order of Young's moudulus random field;

    %   t: Thickness field;
    %   eole_order_t: EOLE order of thickness random field;

    %   nu: Poisson's ratio (univariate in the whole domain);

    %   boundary_condition: Structure array for boundary conditions;
    %       eg:
    %            allDOFs = 1:2*(nely + 1)*(nelx + 1);
    %            boundary_condition.fixedDOFs = 1:1:2*(nely + 1); % Cantilever beam
    %            boundary_condition.freeDOFs = setdiff(allDOFs, boundary_condition.fixedDOFs);

    %   f_components:
    %       f_components.n: Number of loads;
    %       f_components.DOF: DOFs of loads;
    %       f_components.nominal_value: Nominal values of loads;

    %       eg:
    %       f_components.n = 2;
    %       f_components.DOF = [1; 2; 3; 4]; % Column vector
    %       f_components.nominal_value = [10; -10; 20; -20]; % Column vector
    %
    %       f = sparse(f_components.DOF, ones(2*f_components.n, 1), f_components.nominal_value, 2*(nelx + 1)*(nely + 1), 1);
    %   rand_numbr_f: Num. of uncertain loads;

    %   rho: Element density;
    %   penalty: Penalty factor;

    %   w: Eigenvalue and eigenvector decomposition of characteristic matrix w

    %   phi_e: Mode functions of the Young's modulus random field
    %   error_magnitude_e: Error magnitude of Young's modulus.
    %   uniform_dist_boundary_e: Boundaries of the uniform distribution for Young's modulus randomness.

    %   phi_t: Mode functions of the thickness random field
    %   error_magnitude_t: Error magnitude of thickness.
    %   uniform_dist_boundary_t: Boundaries of the uniform distribution for thickness randomness.

    %   xi_all: Aggregate of Gaussian random variables;

    % Output:

    %   u: Displacement
    %   c: Compliance
    %   d_c: Sensitivity of compliance with respect to element density rho (A nely-row by nelx-column matrix)

    %% Mesh information setup
    nelx = mesh_info.nelx;
    nely = mesh_info.nely;
    ik = mesh_info.ik;
    jk = mesh_info.jk;
    element_node_num = mesh_info.element_node_num;

    % Components of element stiffness
    k_components = zeros(64 * nelx * nely, 1);

    % Initialization of element volume
    % v_element = zeros(nelx*nely, 1);
    % d_v_element_temp = zeros(nelx*nely, 1);

    % Initialization of sensitivity of volume
    % d_v = zeros(nely, nelx);

    if nargout > 1
        % Components of element stiffness sensitivity with respect to rho
        d_k_components = zeros(64 * nelx * nely, 1);
    end

    %% Transform matrix form rho to vector form rho
    rho_vector = reshape(rho, nelx * nely, 1);

    %% Creating loads sparse vector
    if rand_numbr_f > 0

        xi_f = xi_all(1:rand_numbr_f);

        % Considering loads uncertainties
        %   f_components:
        %       f_components.n: Number of uncertain loads;
        %       f_components.DOF: DOFs of loads;
        %       f_components.nominal_value: Nominal values of loads;
        [~, ~, f] = f_total(xi_f, f_components, w, nelx, nely);
    else
        % Creating deterministic loads sparse vector
        f = sparse(f_components.DOF, ones(size(f_components.DOF)), f_components.nominal_value, 2 * (nelx + 1) * (nely + 1), 1);
    end

    %% Constructing stiffness matrix (new)

    if eole_order_e > 0

        xi_e = xi_all(rand_numbr_f + 1:rand_numbr_f + eole_order_e);

        % Calculating random Young's modulus field
        [~, delta_e] = random_field_fun(xi_e, phi_e, mesh_info, uniform_dist_boundary_e);
        e_local = e + error_magnitude_e * delta_e; % Young's modulus within this subroutine
    else
        e_local = e;
    end

    if eole_order_t > 0

        xi_t = xi_all(rand_numbr_f + eole_order_e + 1:rand_numbr_f + eole_order_e + eole_order_t);

        % Calculating random thickness field
        [~, delta_t] = random_field_fun(xi_t, phi_t, mesh_info, uniform_dist_boundary_t);
        t_local = t + error_magnitude_t * delta_t; % Element thickness within this subroutine
    else
        t_local = t;
    end

    %%
    for ii = 1:nelx * nely
        idx = element_node_num(ii, :);

        E_Ke = e_local(idx);
        t_Ke = t_local(idx);
        nu_Ke = nu;

        if nargout > 1
            [ke_temp, d_ke_temp] = Ke(E_Ke, t_Ke, nu_Ke, rho_vector(ii), penalty);
            k_components((64 * (ii - 1) + 1):(64 * ii)) = ke_temp(:);
            d_k_components((64 * (ii - 1) + 1):(64 * ii)) = d_ke_temp(:);
        else
            ke_temp = Ke(E_Ke, t_Ke, nu_Ke, rho_vector(ii), penalty);
            k_components((64 * (ii - 1) + 1):(64 * ii)) = ke_temp(:);
        end

    end

    %% Solving for displacement u
    k = sparse(ik, jk, k_components);
    k = (k + k.') / 2; % To ensure symmetry of k

    u = zeros(2 * (nely + 1) * (nelx + 1), 1);

    u(boundary_condition.freeDOFs, :) = ...
        k(boundary_condition.freeDOFs, boundary_condition.freeDOFs) \ f(boundary_condition.freeDOFs, :);

    u(boundary_condition.fixedDOFs, :) = 0;

    %% Calculating sensitivity of compliance with respect to element density
    if nargout > 1

        % Calculating compliance
        c = u.' * k * u;

        u_idx = zeros(1, 8);
        d_c = zeros(1, nelx * nely);

        for ii = 1:nelx * nely
            idx = element_node_num(ii, :);
            u_idx(1:2:7) = 2 * idx - 1;
            u_idx(2:2:8) = 2 * idx;
            dk = d_k_components((64 * (ii - 1) + 1):(64 * ii));
            dk = reshape(dk, 8, 8);
            d_c(ii) =- u(u_idx).' * dk * u(u_idx);
        end

        d_c = reshape(d_c, nely, nelx);

    end

end % End of FEA()
