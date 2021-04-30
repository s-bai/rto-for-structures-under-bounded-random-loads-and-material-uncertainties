function pce_gamma = pce_normal_factor(multi_index)
%% This function calculates the normalization factors for PCE
%   Written by Song Bai.

%   Input:
%       multi_index: Multi-index tuple (One dimensional row or column vector)

%   Output:
%       pce_gamma: Normalization factors for PCE 

%   Note: 
%       The built-in Hermite polynomials function in MATLAB is the
%       so-called "Physicists' Hermite Polynomial H_n(x)". Denoting the
%       "Probabilists' Hermite Polynomial" by He_n(x) and one has:
%       He_n(x) = 2^(-n/2)*H_n(x/sqrt(2)).
%       The normal factors herein are calculated using Physicists' Hermite Polynomial
%       by means of variable replacement.

%% The first 10 orders of univariate gamma_i
gamma_i_univariate = ...
    [
    1       % 0-th
    1       % 1-st
    2       % 2-nd
    6       % 3-rd
    24      % 4-th
    120     % 5-th
    720     % 6-th
    5040    % 7-th
    40320   % 8-th
    362880  % 9-th
    3628800 % 10-th
    ];

% Calculating multi-index normal factor
pce_gamma = prod(gamma_i_univariate(multi_index + 1));

end % End of pce_normal_factor()
