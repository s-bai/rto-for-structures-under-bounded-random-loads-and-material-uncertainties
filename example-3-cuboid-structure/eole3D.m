function [eigen_value, phi3D] = eole3D(mesh_info3D, n_eole, correlation_length_factor)
    %%

    %%

    %%
    nelx = mesh_info3D.nelx;
    nely = mesh_info3D.nely;
    nelz = mesh_info3D.nelz;
    xyz = mesh_info3D.xyz;

    %%
    filename = strcat('EOLE Phi 3D_', num2str(nelx), 'by', num2str(nely), 'by', num2str(nelz), '_');
    filename = strcat(filename, 'n_EOLE=', num2str(n_eole), '_lc=', num2str(correlation_length_factor), '.mat');

    if isfile(filename)
        %%
        fprintf('\nExisting Phi data file found!\n');
        data_from_file = load(filename);
        phi3D = data_from_file.phi3D;
        eigen_value = data_from_file.eigen_value;
    else
        %%
        correlation_matrix = zeros(nelx * nely * nelz, nelx * nely * nelz);
        correlation_length = correlation_length_factor * max([nelx, nely, nelz]);

        %%
        for ii = 1:nelx * nely * nelz
            xyz_temp = repmat(xyz(ii, :), nelx * nely * nelz, 1);
            correlation_matrix(ii, :) = exp(-vecnorm((xyz_temp - xyz).').^2 / correlation_length^2);
        end

        %%
        [eigen_vector, eigen_value] = eigs(correlation_matrix, n_eole);
        eigen_value = diag(eigen_value);

        %%
        phi3D = zeros(n_eole, nelx * nely * nelz);

        for ii = 1:n_eole
            phi3D(ii, :) = 1 / sqrt(eigen_value(ii)) * eigen_vector(:, ii).' * correlation_matrix;
        end

        save(filename, 'phi3D', 'eigen_value');
    end

end % End of eole3D()
