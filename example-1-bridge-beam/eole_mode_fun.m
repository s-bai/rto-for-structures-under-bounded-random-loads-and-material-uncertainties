function [phi] = eole_mode_fun(mesh_info, correlation_length_factor, n_eole_order)
    %% Calculating mode functions of the random field
    %
    % Input:
    %   mesh_info: Mesh information;
    %   correlation_length_factor: Correlation length factor;
    %   n_eole_order: Num. of terms truncated for EOLE approximation;
    %
    % Output:
    %   phi: Mode functions (In rows of the matrix);
    %
    % The correlation function is chosen as
    %                           ||x2 - x1||^2
    %       c(x1, x2) = exp( - ---------------- )
    %                                lc^2

    %%
    nelx = mesh_info.nelx;
    nely = mesh_info.nely;
    xy = mesh_info.xy;

    %% Searching for existing mode function file
    nelx_string = num2str(nelx);
    nely_string = num2str(nely);
    CorLenStr = num2str(correlation_length_factor);
    n_eole_order_string = num2str(n_eole_order);
    filename = strjoin({'Mode functions of mesh', nelx_string, 'by', nely_string, ...
                        'CorreLen factor =', CorLenStr, ...
                        'EOLE terms =', n_eole_order_string});
    filename = join([filename, '.mat']);
    searchresult = dir(filename);

    if size(searchresult, 1) > 0
        %% Load file
        disp('Existing MAT File found:');
        disp(strcat('"', filename, '"'));
        load(filename, 'phi');
        fprintf('\n');
    else
        %% Calculating correlation matrix
        % Cross correlation matrix of the i-th sample point is the corresponging
        % i-th column of correlation matrix.
        correlation_matrix = zeros((nelx + 1) * (nely + 1), (nelx + 1) * (nely + 1));
        correlation_length = correlation_length_factor * max(nelx, nely);

        %%
        for ii = 1:(nelx + 1) * (nely + 1)
            xy_temp = repmat(xy(ii, :), (nelx + 1) * (nely + 1), 1);
            correlation_matrix(ii, :) = exp(-vecnorm((xy_temp - xy).').^2 / correlation_length^2);
        end

        %% Calculate eigenvectors and eigenvalues of the correlation matrix
        % Calculating the first n_eole_order eigenvalues/vectors;
        % eigenV: Eigenvectors;
        % eigenD: Eigenvalues;
        [eigenV, eigenD] = eigs(correlation_matrix, n_eole_order);
        eigenD = diag(eigenD); % Put eigenvalues into column vector;

        %% Discretization of the underlying Gaussian random field using EOLE method
        % Calculating mode functions (Stored in rows of the obtained matrix)
        phi = zeros(n_eole_order, (nelx + 1) * (nely + 1));

        for ii = 1:n_eole_order
            phi(ii, :) = 1 / sqrt(eigenD(ii)) * eigenV(:, ii).' * correlation_matrix;
        end

        %% Save file
        save(filename, 'phi');
    end

end % End of function eole_mode_fun()
