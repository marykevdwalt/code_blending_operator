function [f_y_j] = function_evaluate(y, f, j)

% given an array y = [y(1),...,y(n)]
% given an array of function values f = [f(1),...,f(L)] evaluated on linspace(y(1), y(n), L) (L very large)
% determine f's value at y(j): f_y_j = f(y(j))

n = length(y);
L = length(f);
dy = linspace(y(1), y(end), L);

% increment is the interval between y-values where f is evaluated
%incr = (y(end) - y(1))/(L-1);
incr = dy(2)-dy(1);

f_y_j = f(round( 1 + ( y(j) - y(1) ) / incr));

end