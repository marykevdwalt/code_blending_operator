function [N_xmk] = b_spline_evaluate_wo_boundaries(x_o, m, L, k)

% x_o is the original knot sequence, determining the interval on which the spline is evaluated
% k can be between 0 and N-m

% dx is the x-values where we will evaluate the b-spline
%dx = linspace(x_o(alpha+1), x_o(beta+2), L);
dx = linspace(x_o(1), x_o(end), L);

% we create a vector containing values of N_{m,x,k}(t)(0<=k<=N-m)
% each row evaluates the b-splines at a different x-value between x_o(alpha+1) and x_o(beta+2)
N_xmk = zeros(L, 1);

for l = 1: L
    N_xmk(l) = b_spline_wo_boundaries(x_o, m, k, dx(l)) ;
end

end