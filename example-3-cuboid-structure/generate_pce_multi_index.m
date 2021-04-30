function [dim_pce, pce_multi_index] = generate_pce_multi_index(n_random_var, n_pce_order)
    %% Generating multi-index for Polynomial Chaos Expansion
    %   Written by Song Bai.
    %   This code is based on the code 'nsumk' by Peter Cotton
    %   https://www.mathworks.com/matlabcentral/fileexchange/28340-nsumk
    %
    %   Input:
    %       n_var: Num. of random variables;
    %       n_pce_degree: Degree of PCE.
    %   Output:
    %       dim_pce: Dimension of PCE;
    %       pce_multi_index: Multi-index of PCE.

    dim_pce = nchoosek(n_random_var + n_pce_order, n_pce_order);
    pce_multi_index = zeros(dim_pce, n_random_var);
    row_index = 2;

    if n_random_var == 1
        pce_multi_index = (0:(dim_pce - 1))';
    else

        for ii = 1:n_pce_order

            [temp_counter, temp_tuple] = nsumk(n_random_var, ii);
            temp_tuple = flipud(temp_tuple);
            pce_multi_index(row_index:(row_index + temp_counter - 1), :) = temp_tuple;
            row_index = row_index + temp_counter;

        end

    end

end % End of generate_pce_multi_index()
