%% Batch submission parameters

nelx = 200;
nely = 100;
volfrac = 0.3;
penalty = 3;
rmin = 2;
k_level = 6;
n_tp_grid = 3;
n_MC_samples = 0;

%%
for kk = 1:3

    rto_2d_batch(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid, n_MC_samples, kk)

end
