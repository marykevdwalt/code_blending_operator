function [J] = interpol_molecule_odd(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L)

% we create a matrix J that contains the molecules J_{i}(x) = sum_{k=0}^{2m-2} b_{i,k} N_{m,t,i+k-(m-1)}(x), i = 0,...,N
% each column corresponds to i=0,...,N
% each row is an evaluation of J_{i}(x) at a different x-value between x_ext(1)=t(1,i+1) and x_ext(end)=t(end,i+1)

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

J = zeros(L, N+1);

t = fine_knot_seq_odd(x, x_o, x_ext, m, N);
b = interpol_molecule_coeff_odd(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L);


%for i = 0:N
for i=alpha-m+1:beta+m-1

	t_i = t(:, i+1)';	
	N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
	
	J(:,i+1) = N_mtki(:, i-(m-1)+m : i+m-1+m) * b(:,i+1);

end


  dx = linspace(x_o(alpha+1), x_o(beta+2), L);

  % plot(dx, N_mtki(:,4), 'g', 'LineWidth', 1);
  % hold on;
  % plot(dx, N_mtki(:,5), 'y', 'LineWidth', 1);
  % hold on;
  % plot(dx, N_mtki(:,6), 'k', 'LineWidth', 1);
  % hold on;
  % plot(dx, N_mtki(:,7), 'r', 'LineWidth', 1);
  % hold on;
  % plot(dx, N_mtki(:,8), 'b', 'LineWidth', 1);

 % plot(dx, J(:, 3), 'k', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 4), 'r', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 5), 'b', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 6), 'g', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 7), 'm', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 8), 'y', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 9), 'k', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 10), 'c', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 11), 'r', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 12), 'b', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 13), 'g', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 14), 'm', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 15), 'y', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 16), 'k', 'LineWidth', 1);
 % hold on;
 % plot(dx, J(:, 17), 'c', 'LineWidth', 1);
 % hold on;

end