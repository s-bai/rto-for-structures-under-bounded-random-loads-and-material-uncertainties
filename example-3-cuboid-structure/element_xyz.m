function mesh_info3D = element_xyz(nelx, nely, nelz)
    %% Meshing for FEA
    % Calculating element coordinates
    % Based on the paper by Liu and Tovar

    %%
    x_coordinate = 0:(nelx - 1);
    y_coordinate = fliplr(0:(nely - 1));
    z_coordinate = 0:(nelz - 1);

    [x, y, z] = meshgrid(x_coordinate, y_coordinate, z_coordinate);

    mesh_info3D.xyz = [x(:), y(:), z(:)];
    mesh_info3D.nelx = nelx;
    mesh_info3D.nely = nely;
    mesh_info3D.nelz = nelz;

    % for ii = 1:nelx*nely*nelz
    %     fprintf('Element #%d: %d %d %d\n', ii, mesh_info3D.xyz(ii, 1), mesh_info3D.xyz(ii, 2), mesh_info3D.xyz(ii, 3));
    % end

end % End of element_xyz()
