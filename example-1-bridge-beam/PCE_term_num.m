function [P] = PCE_term_num(KL_order_M, PCE_order_P)
% This subroutine calculates the number of terms P of PCE
% KL_order_M: Number of terms of K-L expansion
% PCE_order_P: Order of PCE
P_iter = zeros(PCE_order_P, 1);
for ii = 1:PCE_order_P 
P_iter(ii, 1) = factorial(KL_order_M + ii - 1)/...
    (factorial(ii)*factorial(KL_order_M - 1));
end
P = sum(P_iter);
