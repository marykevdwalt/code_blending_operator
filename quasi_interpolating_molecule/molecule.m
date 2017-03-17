function [M] = molecule(y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create a matrix M that contains the molecules M_{i}(x) = sum_{j=0}^{m-1} a_{i,j} N_{m,x,i+j-(m-1)}(x)
% each column corresponds to i
% each row is an evaluation of M_{i}(x) at a different x-value between x_o{alpha} and x_o{beta+1}

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

M = zeros(L, N+2*m-1);

[a] = molecule_coeff(y_ext, x_ext, alpha, beta, m, N, knots);
N_mxk = b_spline_evaluate(x_o, x_ext, alpha, beta, m, N, L);

if strcmp(knots, 'equal')
	% interior molecules
	for i=0:N
		M(:,i+m) = N_mxk(:,i-(m-1)+m : i+m) * a(:,i+m);
	end
	% boundary molecules
    for ll=1:m-1
        M(:,m-ll) = N_mxk(:,-(m-1)+m : -ll+m) * a(0+1:m-1-ll+1, m-ll);
    end
    for r=0:m-2
        M(:,N+1+r+m) = N_mxk(:,N+r-m+2+m : N+m) * a(r+1:m-2+1, N+1+r+m);
    end
	
% 	dx = linspace(x_o(alpha+1), x_o(beta+2), L);
% 	figure
%   plot(dx, (M(:, 1)), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 2)), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 3)), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 4)), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 5)), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 6)), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 7)), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 8)), 'k', 'LineWidth', 1);
% 	hold on;
%   plot(dx, (M(:, end-6)), 'y', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-5)), 'k', 'LineWidth', 1);
%   hold on;
% 	plot(dx, (M(:, end-4)), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-3)), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-2)), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-1)), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end)), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(x_o, zeros(size(x_o)), '*', 'markersize', 3); hold on;
% 	grid minor;
% 	title('Quasi-interpolating molecules');

end

if strcmp(knots, 'half')
	% interior molecules
	for i=0:N
		M(:,i+m) = N_mxk(:,i-(m-1)+m : i+m) * a(:,i+m);
	end
	% boundary molecules
    for ll=1:m-1
        M(:,m-ll) = N_mxk(:,-(m-1)+m : -ll+m) * a(0+1:m-1-ll+1, m-ll);
    end
    for r=1:m-1
        M(:,N+r+m) = N_mxk(:,N+r-m+1+m : N+m) * a(r+1:m-1+1, N+r+m);
    end
	
% 	dx = linspace(y_o(alpha+1), y_o(beta+2), L);
% 	figure
%   plot(dx, (M(:, 1)), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 2)), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 3)), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 4)), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 5)), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 6)), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 7)), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, 8)), 'k', 'LineWidth', 1);
% 	hold on;
%   plot(dx, (M(:, end-10)), 'r', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-9)), 'm', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-8)), 'c', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-7)), 'g', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-6)), 'y', 'LineWidth', 1);
%   hold on;
%   plot(dx, (M(:, end-5)), 'k', 'LineWidth', 1);
%   hold on;
% 	plot(dx, (M(:, end-4)), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-3)), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-2)), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end-1)), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, (M(:, end)), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(x_o, zeros(size(x_o)), '*', 'markersize', 3); hold on;
%   plot(y_o, zeros(size(y_o)), 'o', 'markersize', 3); hold on;
% 	grid minor;
% 	title('Quasi-interpolating molecules');

end

end