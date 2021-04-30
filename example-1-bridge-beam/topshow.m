function topshow(nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid, ft)
%% This function shows the evolution of topology layout
% Input:
%   nelx, nely, volfrac, penalty: For loading result file
%   Deterministic design: ft = 0
%   Robust design: ft = 1

%% Load result file
switch ft
    case 0
        cmd_name = 'det_rto_against_hybrid_uncertainty';
        cmd_parameters = [nelx, nely, volfrac, penalty, rmin];
        iteration_file_name = cmd_history(cmd_name, cmd_parameters);
        det_iter_history = load(iteration_file_name);
        
        max_loop = det_iter_history.loop;
        
    case 1
        cmd_name = 'rto_against_hybrid_uncertainty';
        cmd_parameters = [nelx, nely, volfrac, penalty, rmin, k_level, n_tp_grid];
        iteration_file_name = cmd_history(cmd_name, cmd_parameters);
        rto_iter_history = load(iteration_file_name);
        
        try
            max_loop = rto_iter_history.iter;
        catch
            temp = ...
                find(rto_iter_history.worst_case_compliance_history(:, 1));
            max_loop = temp(end);
        end
end



% try
%     n_ellipse = rdo_iter_history.n_ellipse;
% catch
%     prompt = '# of ellipse model is missing, manually input:';
%     n_ellipse = input(prompt);
% end

% nonzero_iter_step = find(rho_history);
% max_iter_step = nonzero_iter_step(end);

ss = 1;
% rho = zeros();

figure;

switch ft
    case 0
        while ss <= det_iter_history.loop + 1
            
            
            rho_optimal = det_iter_history.x_history(:, :, ss);
            
            
            colormap(gray);
            imagesc(flipud(1 - rho_optimal));
            caxis([0 1]);
            axis equal;
            axis off;
            drawnow;
            
%             pause(0.1);
            ss = ss + 1;
        end
        
    case 1
        while ss <= max_loop + 1
            
            
            rho_iter = rto_iter_history.rho_history(:, :, ss);
            
            colormap(gray);
            imagesc(flipud(1 - rho_iter));
            caxis([0 1]);
            axis equal;
            axis off;
            drawnow;
            
%             pause(0.1);
            ss = ss + 1;
        end
        
        
end















