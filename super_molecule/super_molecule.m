function [Mstar] = super_molecule(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create a matrix Mstar that contains the super-molecules M*_{i}(x) = M_{i}(x) + J_{i}(x) - sum_j M_{i}(y_j) J_{j}(x)
% each column corresponds to i
% each row is an evaluation of M*_{i}(x) at a different x-value between x_o{alpha} and x_o{beta+1}

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% f_o = original function values at sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

M = molecule(y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots);

Mstar = zeros(L, N+2*m-1);

if strcmp(knots, 'equal')
    J = interpol_molecule_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots);
    K = cancel_molecule(J, y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots);
	for i=-m+1:N+m-1
		Mstar(:,i+m) = M(:,i+m) + J(:,i+m) - K(:,i+m);
	end
% 	dx = linspace(x_o(1), x_o(end), L);
% 	figure
% 	plot(dx, Mstar(:, 1), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 2), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 3), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 4), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 5), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 6), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 7), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 8), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 9), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 10), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 11), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 12), 'g', 'LineWidth', 1);
% 	hold on;
%     plot(dx, Mstar(:, 13), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 14), 'k', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 15), 'm', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 16), 'y', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 17), 'c', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 18), 'b', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 19), 'g', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 20), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 21), 'k', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 22), 'c', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 23), 'y', 'LineWidth', 1); %N=19
%     hold on;
% %    plot(dx, Mstar(:, end-6), 'k', 'LineWidth', 1);
% %    hold on;
%     plot(dx, Mstar(:, end-5), 'y', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-4), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-3), 'b', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-2), 'm', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end), 'k', 'LineWidth', 1);
%     hold on;
% 	plot(x_o, ones(size(x_o)), '*', 'markersize', 3);
% 	hold on;
% 	title('Blending i-molecules');
% 	grid minor;

end

if strcmp(knots, 'half')
    J = interpol_molecule_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots);
    K = cancel_molecule(J, y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots);
	for i=-m+1:N+m-1
		Mstar(:,i+m) = M(:,i+m) + J(:,i+m) - K(:,i+m);
	end
% 	dx = linspace(x_o(1), x_o(end), L);
% 	figure
% 	plot(dx, Mstar(:, 1), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 2), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 3), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 4), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 5), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 6), 'g', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 7), 'r', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 8), 'k', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 9), 'm', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 10), 'c', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 11), 'b', 'LineWidth', 1);
% 	hold on;
% 	plot(dx, Mstar(:, 12), 'g', 'LineWidth', 1);
% 	hold on;
%     plot(dx, Mstar(:, 13), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 14), 'k', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 15), 'm', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 16), 'y', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 17), 'c', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 18), 'b', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 19), 'g', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 20), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 21), 'k', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 22), 'c', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 23), 'y', 'LineWidth', 1); %N=19
%     hold on;
%     plot(dx, Mstar(:, 24), 'm', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 25), 'b', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, 26), 'g', 'LineWidth', 1);
%     hold on;
% %    plot(dx, Mstar(:, end-6), 'k', 'LineWidth', 1);
% %    hold on;
%     plot(dx, Mstar(:, end-5), 'y', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-4), 'r', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-3), 'b', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-2), 'm', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end-1), 'g', 'LineWidth', 1);
%     hold on;
%     plot(dx, Mstar(:, end), 'k', 'LineWidth', 1);
%     hold on;
% 	plot(x_o, ones(size(x_o)), '*', 'markersize', 3);
% 	hold on;
%     plot(y_o, ones(size(y_o)), 'o', 'markersize', 3);
% 	title('Blending i-molecules');
% 	grid minor;

end

end