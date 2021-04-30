function [place_holder, v, Dv] = total_volume(mesh_info, t, rho, ...
        eole_order_t, xi_all, phi_t, error_magnitude_t, uniform_dist_boundary_t)
    %% Calculate total material volume and its sensitivity

    % Input:
    %   mesh_info
    %       mesh_info.nelx: Num. of elements in direction x;
    %       mesh_info.nely: Num. of elements in direction y;
    %       mesh_info.element_node_num: Element node numbers;
    %   t: Nominal thickness field ;
    %   rho: Density matrix;

    %   eole_order_t: Num. of random variables for thickness random field;

    %   xi_all: Total random variable vector;

    % Output:
    %   placeholder: Placeholder for output consistency with FEA subroutine;
    %   v: Total volume
    %   Dv: Sensitivity of volume with respect to element density.

    %% Ininitalization
    nelx = mesh_info.nelx;
    nely = mesh_info.nely;
    element_node_num = mesh_info.element_node_num;
    v_element = zeros(nelx * nely, 1);
    Dv_element_temp = zeros(nelx * nely, 1);

    %% Retrieving thickness random variables from total random variable vector
    xi_t = xi_all(end - eole_order_t + 1:end);

    %% Calculating perturbated thickness field
    [~, delta_t] = random_field_fun(xi_t, phi_t, mesh_info, uniform_dist_boundary_t);
    t_local = t + error_magnitude_t * delta_t; % Element thickness within this subroutine

    %% Transform matrix form rho to vector form rho
    rho_vector = reshape(rho, nelx * nely, 1);

    %%
    for ii = 1:nelx * nely
        idx = element_node_num(ii, :);

        t_Ke = t_local(idx);
        [v_element(ii), Dv_element_temp(ii)] = element_volume(t_Ke, rho_vector(ii));
    end

    v = sum(v_element);
    Dv = reshape(Dv_element_temp, nely, nelx);

    place_holder = 0;

end % End of total_volume()
