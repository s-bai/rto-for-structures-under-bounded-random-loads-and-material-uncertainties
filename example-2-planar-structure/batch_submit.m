%% Batch submission parameters

nelx = 180;
nely = 120;
volfrac = 0.3;
penalty = 3;
rmin = 2.4;
k_level = 5;
n_tp_grid = 3;
n_MC_samples = 0;

%%
for kk = 1:5

    rto_2d_batch(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid, n_MC_samples, kk)

end
