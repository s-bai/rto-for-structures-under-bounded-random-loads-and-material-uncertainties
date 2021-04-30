function rto_loads3D(nelx, nely, nelz, volfrac, penal, rmin, k_level, n_pce_order)
%% Robust topology optimization against uncertainties for 3D structures
%  Written by Song Bai, supervised by Prof. Zhan Kang
%  Based on 'AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR'

%% USER-DEFINED LOOP PARAMETERS
max_iter = 100;    % Maximum number of iterations
min_change = 0.01;      % Terminarion criterion
displayflag = 1;  % Display structure flag
beta0 = 1;
n_rand_var = 3;


%% USER-DEFINED MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio


%% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);




%% USER-DEFINED LOAD DOFs
%% Cantilever beam at right lower edge
% [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
% loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
% loaddof = 3*loadnid(:) - 1;                             % DOFs

%% Cube beam with concentrate loads +Fx at right center
% Coordinates
[il,jl,kl] = meshgrid(nelx, nely/2, nelz/2);
% Node IDs
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);
% DOFs
loaddof = [3*loadnid(:) - 2, 3*loadnid(:) - 1, 3*loadnid(:)];



%%
% Characteristic matrix
f_components.w = [3, 0, 0; 0, 3, 0; 0, 0, 1];
% Eigen decomposition
[f_components.Q, f_components.SIGMA] = eig(f_components.w);
%
f_components.T = f_components.Q*sqrt(f_components.SIGMA)*f_components.Q.';
%
f_components.T_inverse = inv(f_components.T);
% Num. of uncertain loads
f_components.n = 1;
%
% DOFs
f_components.DOF = loaddof;
% Nominal value
f_components.f0 = [1; 0; 0];
% Total DOFs
f_components.total_DOFs = ndof;





%% USER-DEFINED SUPPORT FIXED DOFs
%% Cantilever beam fixed at left facet
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs

freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);


%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);


%% INITIALIZE ITERATION
rho = repmat(volfrac,[nely,nelx,nelz]);
rho_Phys = rho;
iter = 0;
change = 1;

%% Mesh information
mesh_info.nelx = nelx;
mesh_info.nely = nely;
mesh_info.nelz = nelz;
mesh_info.ndof = ndof;
mesh_info.nele = nele;
mesh_info.iK = iK;
mesh_info.jK = jK;
mesh_info.freedofs = freedofs;

%% Material properties
material_properties.edofMat = edofMat;
material_properties.KE = KE;
material_properties.Emin = Emin;
material_properties.E0 = E0*ones(nely, nelx, nelz);


%% Generating Sparse Grid quadrature nodes and the corresponding weights
[SG_nodes, SG_weights] = nwspgr('KPN', n_rand_var, k_level);


%% Generating PCE multi-index
[PCE_dimension, PCE_multi_index] = generate_pce_multi_index(n_rand_var, n_pce_order);


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


%% Iteration histories
iter_history = zeros(max_iter, 3);
rho_history = zeros(nely, nelx, nelz, max_iter + 1);
rho_history(:, :, :, 1) = rho;

%% START ITERATION
while change > min_change && iter < max_iter
    iter = iter+1;
    
    %% Creating function handle for FEA
    FEA_handle = @(xi_all) FEA3D(material_properties, mesh_info, rho_Phys, penal, f_components, xi_all);
    
    
    
    %% Calculating structural response at SPARSE GRID quadrature nodes
    [c_nodal_value_SG, Dc_nodal_value_SG] = f_integration_nodal3D(FEA_handle, SG_nodes, nelx, nely, nelz);
    
    
    %% Calculating PCE coefficients
    for ii = 1:PCE_dimension
        c_hat_SG(ii, 1) = (1/PCE_gamma(ii))*sum(c_nodal_value_SG.*PCE_base_nodal_value_SG(:, ii).*SG_weights);
        
        Dc_nodal_value_SG_temp = Dc_nodal_value_SG;
        for jj = 1:size(SG_nodes, 1)
            Dc_nodal_value_SG_temp(:, :, :, jj) = ...
                Dc_nodal_value_SG_temp(:, :, :, jj)*PCE_base_nodal_value_SG(jj, ii)*SG_weights(jj, 1);
        end
        Dc_hat_SG(:, :, :, ii) = (1/PCE_gamma(ii, 1))*sum(Dc_nodal_value_SG_temp, 4);
    end
    
    
    %% Mean value and the mean value sensitivity of compliance (SPARSE GRID)
    c_mean_PCE_SG = c_hat_SG(1, 1);
    Dc_mean_PCE_SG = Dc_hat_SG(:, :, :, 1);
    
    
    %% STD and STD sensitivity of compliance (SPARSE GRID)
    c_STD_PCE_SG = sqrt(PCE_gamma(2:end, 1).'*(c_hat_SG(2:end).^2));
    
    Dc_STD_PCE_SG_temp = Dc_hat_SG;
    for ii = 1:PCE_dimension
        Dc_STD_PCE_SG_temp(:, :, :, ii) = PCE_gamma(ii, 1)*c_hat_SG(ii, 1)*Dc_hat_SG(:, :, :, ii);
    end
    Dc_STD_PCE_SG = (1/c_STD_PCE_SG)*sum(Dc_STD_PCE_SG_temp, 4);
    
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    obj = c_mean_PCE_SG + beta0*c_STD_PCE_SG;
    Dobj = Dc_mean_PCE_SG + beta0*Dc_STD_PCE_SG;
    
    iter_history(iter, 1:3) = [c_mean_PCE_SG, c_STD_PCE_SG, obj];
    
    
    %% FILTERING AND MODIFICATION OF SENSITIVITIES
    Dv = ones(nely,nelx,nelz);
    Dobj(:) = H*(Dobj(:)./Hs);
    Dv(:) = H*(Dv(:)./Hs);
    
    
    %% OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        rho_new = max(0,max(rho-move,min(1,min(rho+move,rho.*sqrt(-Dobj./Dv/lmid)))));
        rho_Phys(:) = (H*rho_new(:))./Hs;
        if sum(rho_Phys(:)) > volfrac*nele
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
    change = max(abs(rho_new(:)-rho(:)));
    rho = rho_new;
    
    rho_history(:, :, :, iter + 1) = rho_Phys;
    
    
    %% PRINT RESULTS
    fprintf('\nMean   STD   Objective\n');
    disp([c_mean_PCE_SG, c_STD_PCE_SG, obj]);
    fprintf('\n It.:%5i  Vol.:%7.3f  ch.:%7.3f\n', iter, mean(rho_Phys(:)), change);
    
    
    %% PLOT DENSITIES
    if displayflag, clf; display_3D(rho_Phys); end % #ok<UNRCH>
    
    %%
    save('rto3D.mat', 'rho_history', 'iter_history');
    
    
end

fprintf('\nSolution converged after %d iteration steps.\n', iter);

clf; display_3D(rho_Phys);


end % End of rto_loads3D



