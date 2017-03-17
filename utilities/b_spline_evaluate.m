function [N_xmk] = b_spline_evaluate(x_o, x_ext, alpha, beta, m, N, L)

% x_o is the original knot sequence, determining the interval on which the spline is evaluated
% x_ext is the knot sequence to be used
% alpha is the index of the first knot in x_o where we will evaluate the b-splines
% beta+1 is the index of the last knot in x_o where we will evaluate the b-splines

% dx is the x-values where we will evaluate the b-spline
dx = linspace(x_o(1), x_o(end), L);

% we create a matrix containing b-spline values
% each column contains a b-spline starting at a different knot (N_{m,x,-(m-1)},...,N_{m,x,N})
% each row evaluates the b-splines at a different x-value between x_o(alpha+1) and x_o(beta+2)
N_xmk = zeros(L, m+N);

for k = -(m-1):1:N

	for l = 1: L
		N_xmk(l, k+m) = b_spline(x_ext, m, N, k, dx(l)) ;
	end
	
%	N_xmk(:, k+m) = N_xmk(:, k+m) ./ norm(N_xmk(:, k+m)) ;

end

% figure
%  plot(dx, N_xmk(:,1), 'g', 'LineWidth', 1);
%  hold on;
%  plot(dx, N_xmk(:,2), 'y', 'LineWidth', 1);
%  hold on;
%  plot(dx, N_xmk(:,3), 'k', 'LineWidth', 1);
%  hold on;
%  plot(dx, N_xmk(:,4), 'r', 'LineWidth', 1);
%  hold on;
%  plot(dx, N_xmk(:,5), 'm', 'LineWidth', 1);
%  hold on;
%  plot(dx, N_xmk(:,6), 'c', 'LineWidth', 1);
%  hold on;
% 
% plot(dx, N_xmk(:,7), 'g', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,8), 'm', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,9), 'y', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,10), 'b', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,11), 'r', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,12), 'g', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,13), 'm', 'LineWidth', 1);
% hold on;
% 
% plot(dx, N_xmk(:,17), 'b', 'LineWidth', 1);
% hold on
% plot(dx, N_xmk(:,18), 'k', 'LineWidth', 1);
% hold on
% plot(dx, N_xmk(:,19), 'c', 'LineWidth', 1);
% hold on
% plot(dx, N_xmk(:,20), 'm', 'LineWidth', 1);
% hold on
% plot(dx, N_xmk(:,21), 'g', 'LineWidth', 1);
% 
% plot(dx, N_xmk(:,end-4), 'b', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,end-3), 'm', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,end-2), 'k', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,end-1), 'r', 'LineWidth', 1);
% hold on;
% plot(dx, N_xmk(:,end), 'g', 'LineWidth', 1);

end