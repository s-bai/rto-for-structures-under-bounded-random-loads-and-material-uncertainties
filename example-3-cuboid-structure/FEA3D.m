function [U, c, Dc] = ...
        FEA3D(material_properties, mesh_info3D, rho_Phys, penal, ...
        f_components, n_rand_f, eole_order_e, err_magnitude_e, phi3D, xi_all)
    %% FEA subroutine for three-dimensional structures
    % Based on the code by Liu and Tovar
    % Written by Song Bai, supervised by Prof. Zhan Kamg

    %% Input:
    %   rho_Phys: Element physical densities;
    %   f_components: Components of loads;
    %   xi_all: Total random variable row vector;
    %   phi3D: Elemental EOLE mode functions

    %% Output:
    %   U: Displacement;
    %   c: Mean compliance;
    %   Dc: Sensitivities of mean compliance with respect to element densities;

    %% Mesh information
    nelx = mesh_info3D.nelx;
    nely = mesh_info3D.nely;
    nelz = mesh_info3D.nelz;
    ndof = mesh_info3D.ndof;
    nele = mesh_info3D.nele;
    iK = mesh_info3D.iK;
    jK = mesh_info3D.jK;
    freedofs = mesh_info3D.freedofs;

    %% Material properties
    edofMat = material_properties.edofMat;
    KE = material_properties.KE;
    Emin = material_properties.Emin;
    E0 = material_properties.E0;

    %%
    if n_rand_f > 0
        xi_all_f = xi_all(1:n_rand_f);
        [~, f] = f_total3d(xi_all_f, f_components);
    else
        f = f_components.f0;
    end

    if eole_order_e > 0
        % size(x_all_e) = (1, n_eole_e)
        xi_all_e = xi_all((n_rand_f + 1):end);

        % Calculating random Young's modulus field
        [~, E0_delta] = gaussian_to_uniform(xi_all_e, phi3D);

        E0_column = E0(:) + err_magnitude_e * E0_delta.';
        E0 = reshape(E0_column, [nely, nelx, nelz]);
    end

    %%
    U = zeros(ndof, 1);

    sK = reshape(KE(:) * (Emin + rho_Phys(:).'.^penal .* (E0_column - Emin).'), 24 * 24 * nele, 1);
    % sK = KE(:)*(Emin+rho_Phys(:).'.^penal.*(E0_column - Emin).');

    K = sparse(iK, jK, sK); K = (K + K') / 2;
    
    %% Direct method
    % U(freedofs,:) = K(freedofs,freedofs)\f(freedofs,:);
    
    %% PCG solver
    tolit = 1e-5;
    maxit = 8000;
    M = diag(diag(K(freedofs, freedofs)));
    [U(freedofs, :), ~] = pcg(K(freedofs, freedofs), f(freedofs, :), tolit, maxit, M);

    ce = reshape(sum((U(edofMat) * KE) .* U(edofMat), 2), [nely, nelx, nelz]);
    c = sum(sum(sum((Emin + rho_Phys.^penal .* (E0 - Emin)) .* ce)));
    Dc = -penal * (E0 - Emin) .* rho_Phys.^(penal - 1) .* ce;

end % End of FEA3D()
