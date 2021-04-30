%% Batch submission parameters for rto_3D
nelx = 20;
nely = 40;
nelz = 40;
volfrac = 0.2;
penalty = 3;
rmin = 1.5;
n_tp_grid = 3;
n_MC_samples = 0;

%% Load data from Excel file
% data_table = readtable('rto_3D.xlsx', 'PreserveVariableNames', true);
data_table = readtable('rto_3D.xlsx');
kk = data_table.serial_number(end);

%%
for sn = 1:1

    fprintf('\nJob #%d/%d\n', sn, kk);

    rto_3D(nelx, nely, nelz, volfrac, penalty, rmin, n_tp_grid, n_MC_samples, sn);

end

fprintf('\nAll jobs done.\n');
