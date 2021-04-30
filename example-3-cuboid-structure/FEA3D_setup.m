function [mesh_info, material_properties] = FEA3D_setup(nelx, nely, nelz, il, jl, kl, iif, jf, kf, E0, Emin, nu, w, f0)
    %% Parameters setup for FEA3D

    %% PREPARE FINITE ELEMENT ANALYSIS
    nele = nelx * nely * nelz;
    ndof = 3 * (nelx + 1) * (nely + 1) * (nelz + 1);

    %% USER-DEFINED LOAD DOFs
    % Node IDs
    loadnid = kl * (nelx + 1) * (nely + 1) + il * (nely + 1) + (nely + 1 - jl);
    % DOFs
    loaddof = [3 * loadnid(:) - 2, 3 * loadnid(:) - 1, 3 * loadnid(:)];

    %%
    % Characteristic matrix
    f_components.w = w;
    % Eigen decomposition
    [f_components.Q, f_components.SIGMA] = eig(f_components.w);
    %
    f_components.T = f_components.Q * sqrt(f_components.SIGMA) * f_components.Q.';
    %
    f_components.T_inverse = inv(f_components.T);
    % Num. of uncertain loads
    f_components.n = 1;
    % DOFs
    f_components.DOF = loaddof;
    % Nominal value
    f_components.f0 = f0;
    % Total DOFs
    f_components.total_DOFs = ndof;

    %% USER-DEFINED SUPPORT FIXED DOFs
    fixednid = kf * (nelx + 1) * (nely + 1) + iif * (nely + 1) + (nely + 1 - jf); % Node IDs
    fixeddof = [3 * fixednid(:); 3 * fixednid(:) - 1; 3 * fixednid(:) - 2]; % DOFs
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

    %% Mesh information
    mesh_info.nelx = nelx;
    mesh_info.nelx = nely;
    mesh_info.nelx = nelz;
    mesh_info.ndof = ndof;
    mesh_info.nele = nele;
    mesh_info.iK = iK;
    mesh_info.jK = jK;
    mesh_info.freedofs = freedofs;

    %% Material properties
    material_properties.edofMat = edofMat;
    material_properties.KE = KE;
    material_properties.Emin = Emin;
    material_properties.E0 = E0;

end % End of FEA3D_setup()
