function [place_holder1, pce_base_value, place_holder2] = pce_base_fun(pce_index, xi_all)
%% PCE base function

%% Input:
%   pce_index:      PCE index;
%   xi_all: Total random variable vector;


%% Output:
%   place_holder1: For output consistency with FEA subroutine;
%   fun_value: PCE base value;
%   place_holder2: For output consistency with FEA subroutine;



%% Note that the Hermite polynomial in MATLAB is the 'physicists' version.
% Variable replacement
xi_all = xi_all./sqrt(2);



hermiteH10 = @(xi) ...
    [1; ...
    2*xi; ...
    4*xi^2 - 2; ...
    8*xi^3 - 12*xi; ...
    16*xi^4 - 48*xi^2 + 12; ...
    32*xi^5 - 160*xi^3 + 120*xi; ...
    64*xi^6 - 480*xi^4 + 720*xi^2 - 120; ...
    128*xi^7 - 1344*xi^5 + 3360*xi^3 - 1680*xi; ...
    256*xi^8 - 3584*xi^6 + 13440*xi^4 - 13440*xi^2 + 1680; ...
    512*xi^9 - 9216*xi^7 + 48384*xi^5 - 80640*xi^3 + 30240*xi; ...
    1024*xi^10 - 23040*xi^8 + 161280*xi^6 - 403200*xi^4 + 302400*xi^2 - 30240];

%%
hermiteH10_values = zeros(11, size(xi_all, 2));
hermite_transformation_power = zeros(size(xi_all, 2), 1);

for ii = 1:size(xi_all, 2)
    hermiteH10_values(:, ii) = hermiteH10(xi_all(ii));
    hermite_transformation_power(ii, 1) = 2^(-pce_index(ii)/2);
end


%%
row = pce_index + 1;
column = 1:size(xi_all, 2);
index = sub2ind(size(hermiteH10_values), row, column);


%% 
pce_base_value = prod(hermite_transformation_power)*prod(hermiteH10_values(index));




place_holder1 = 0;
place_holder2 = 0;

end % End of pce_base_fun()

