eole_info = load('EOLE Phi 3D_20by40by40_n_EOLE=8_lc=1.8.mat');
nelx = 20;
nely = 40;
nelz = 40;
% eole_mode = zeros(nely, nelx, nelz, 8);

%%
filename_prefix = 'mode';
filename_suffix = '.mat';

%%
for ii = 1:8
    
    eole_mode_ii = reshape(eole_info.phi3D(ii, :), [nely, nelx, nelz]);
    
    filename = strcat(filename_prefix, num2str(ii), filename_suffix);
    save(filename, 'eole_mode_ii');
end