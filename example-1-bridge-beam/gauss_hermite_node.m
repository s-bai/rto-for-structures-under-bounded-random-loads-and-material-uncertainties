function [gh_nodes, gh_weights] = gauss_hermite_node(n)
%% Finding the quadrature nodes for one dimensional Gauss-Hermite quadrature
%   Input:
%       n: Num. of quadrature nodes;
%          n is also the order of Physicists' Hermite polynomial 
%          (n = 3 by default);
%   Output:
%       gh_nodes: Gauss-Hermite quadrature nodes

if nargin < 1
    n = 3;
end

if n <= 0
    fprintf('%s\n', "n MUST be a positive integer!");
    fprintf('%s\n', "n is set to be n = 3 by default.");
    n = 3;
end

syms x;
eqn = hermiteH(n, x) == 0;
solution_sym = vpasolve(eqn, x);
gh_nodes = double(solution_sym);

gh_weights = zeros(n, 1);
for ii = 1:n
    gh_weights(ii) = ( 2^(n-1)*factorial(n)*sqrt(pi) / ...
        (n^2*hermiteH(n - 1, gh_nodes(ii))^2) );
end

end