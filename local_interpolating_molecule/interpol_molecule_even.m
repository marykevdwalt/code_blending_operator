function [J] = interpol_molecule_even(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives)

% we create a matrix J that contains the molecules J_{i}(x) = sum_{k=0}^{2m-3} b_{i,k} N_{m,t,i+k-(m-1)}(x), i = 0,...,N+1 plus derivatices left and right
% each column corresponds to i=0,...,N+1 plus derivatives
% each row is an evaluation of J_{i}(x) at a different x-value between x_ext(1)=t(1,i+1) and x_ext(end)=t(end,i+1)

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

J = zeros(L, N+7);

t = fine_knot_seq_even(x, x_o, x_ext, m, N, boundaries, derivatives);
b = interpol_molecule_coeff_even(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives);

if (strcmp(boundaries, 'left')) && (strcmp(derivatives, 'one'))
	% interior molecules
	for i=m-1:beta+m-1
		t_i = t(:, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:, i-(m-1)+m : i+m-2+m) * b(1:6,i+4);
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		t_i = t(1:end-1, i+3)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+3) = N_mtki(:,1:4) * b(2:5,i+3);
	end
	for i=0:2
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,i-(m-1)+m : i+m-3+m) * b(1:5,i+4);
	end

	% dx = linspace(x_o(1), x_o(beta+2), L);
end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one'))
	% interior molecules
	for i=m-1:N-3
		t_i = t(:, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:, i-(m-1)+m : i+m-2+m) * b(1:6,i+4);
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		t_i = t(1:end-1, i+3)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+3) = N_mtki(:,1:4) * b(2:5,i+3);
	end
	for i=0:2
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,i-(m-1)+m : i+m-3+m) * b(1:5,i+4);
	end
	dx = linspace(x_o(alpha+1), x_o(beta+2), L);
	figure
	plot(dx, J(:,1), 'r'); hold on;
	plot(dx, J(:,2), 'g'); hold on;
	plot(dx, J(:,3), 'b'); hold on;
	plot(dx, J(:,4), 'm'); hold on;
	plot(dx, J(:,5), 'k'); hold on;
	plot(dx, J(:,6), 'c'); hold on;
	plot(dx, J(:,7), 'y'); hold on;
	plot(dx, J(:,8), 'r'); hold on;
	legend('J^3_0', 'J^2_0', 'J^1_0', 'J_0', 'J_1', 'J_2', 'J_3', 'J_4');
	grid minor;

	% boundary molecules (RIGHT)
	for i=N-2
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,end-6 : end-2) * b(1:5,i+4);
	end
	for i=N-1
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,end-5 : end-1) * b(1:5,i+4);
	end
	for i=N
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,end-4 : end) * b(1:5,i+4);
	end
	for i=N+1
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,end-3 : end) * b(1:4,i+4);
	end
	for i=N+1 % first derivative
		t_i = t(1:end-2, i+5)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+5) = N_mtki(:,end-2 : end) * b(1:3,i+5);
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one_zero'))
	% interior molecules
	for i=m-1:N-3
		t_i = t(:, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:, i-(m-1)+m : i+m-2+m) * b(1:6,i+4);
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		t_i = t(1:end-1, i+3)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+3) = N_mtki(:,1:4) * b(2:5,i+3);
	end
	for i=0:2
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,i-(m-1)+m : i+m-3+m) * b(1:5,i+4);
	end

	% boundary molecules (RIGHT)
	for i=N-2
		t_i = t(1:end-2, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+1, L);
		J(:,i+4) = N_mtki(:,end-6 : end-2) * b(1:5,i+4);
	end
	for i=N-1
		t_i = t(1:end-3, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N, L);
		J(:,i+4) = N_mtki(:,end-4 : end-1) * b(1:4,i+4);
	end
	for i=N
		t_i = t(1:end-3, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N, L);
		J(:,i+4) = N_mtki(:,end-3 : end) * b(1:4,i+4);
	end
	for i=N+1
		t_i = t(1:end-3, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N, L);
		J(:,i+4) = N_mtki(:,end-2 : end) * b(1:3,i+4);
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'three_two'))
	% interior molecules
	for i=m:N-3
		t_i = t(1:end-1, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:, i-(m-1)+m : i+m-2+m) * b(1:6,i+4);
	end
	
	% boundary molecules (LEFT)
	% for i=0 %third derivative
		% t_i = t(1:end, 1)';
		% N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+3, L);
		% J(:,1) = N_mtki(:,1:4) * b(1:4,1);
	% end
	% for i=0 %second derivative
		% t_i = t(1:end, 2)';
		% N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+3, L);
		% J(:,2) = N_mtki(:,1:5) * b(1:5,2);
	% end
	% for i=0 %first derivative
		% t_i = t(1:end, 3)';
		% N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+3, L);
		% J(:,3) = N_mtki(:,1:6) * b(1:6,3);
	% end
	[bl,tl,br,tr] = simpler_interpol_molecule_coeff(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives);
	xx = [ x_o(1) + [-m+1:-1]*1e-3 , x , x_o(end) + [1:m-1]*1e-3 ];
	N_mtki = b_spline_evaluate(x_o, xx, alpha, beta, m, 2*N+1, L);
	NN = b_spline_evaluate(x_o, tl{1}, alpha, beta, m, 4, L); N_t00 = NN(:,0 +4);
	NN = b_spline_evaluate(x_o, tl{2}, alpha, beta, m, 4, L); N_t01 = NN(:,-1 +4);
	NN = b_spline_evaluate(x_o, tl{3}, alpha, beta, m, 4, L); N_t02 = NN(:,-2 +4);
	NN = b_spline_evaluate(x_o, tl{4}, alpha, beta, m, 4, L); N_t03 = NN(:,-3 +4);
	N_bdry_left = [N_t00, N_t01, N_t02, N_t03];
	NN = b_spline_evaluate(x_o, tr{1}, alpha, beta, m, length(tr{1})-8, L); N_tNp10 = NN(:,end);
	NN = b_spline_evaluate(x_o, tr{2}, alpha, beta, m, length(tr{2})-8, L); N_tNp11 = NN(:,end-1);
	NN = b_spline_evaluate(x_o, tr{3}, alpha, beta, m, length(tr{3})-8, L); N_tNp12 = NN(:,end-2);
	N_bdry_right = [N_tNp10, N_tNp11, N_tNp12];
	for i=-3:0
		J(:,i+4) = N_bdry_left * bl(:,i+4);
	end

%	for i=0:3
	for i=1:3
		t_i = t(1:end, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+3, L);
		J(:,i+4) = N_mtki(:,i-(m-1)+m : i+m-1+m) * b(1:7,i+4);
	end

	dx = linspace(x_o(alpha+1), x_o(beta+2), L);
	% figure
	% plot(dx, J(:,1), 'r'); hold on;
	% plot(dx, J(:,2), 'g'); hold on;
	% plot(dx, J(:,3), 'b'); hold on;
	% plot(dx, J(:,4), 'm'); hold on;
	% plot(dx, J(:,5), 'k'); hold on;
	% plot(dx, J(:,6), 'c'); hold on;
	% plot(dx, J(:,7), 'y'); hold on;
	% plot(dx, J(:,8), 'r'); hold on;
	% legend('J^3_0', 'J^2_0', 'J^1_0', 'J_0', 'J_1', 'J_2', 'J_3', 'J_4');
	% grid minor;
	
	% boundary molecules (RIGHT)
	for i=N-2:N
		t_i = t(1:end-1, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:,i-(m-1)+m : i+m-2+m) * b(1:6,i+4);
	end
	for i=N+1
		t_i = t(1:end-1, i+4)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,i+4) = N_mtki(:,end-4 : end) * b(1:5,i+4);
	end
	for i=N+1 % first derivative
		t_i = t(1:end-1, N+6)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,N+6) = N_mtki(:,end-3 : end) * b(1:4,N+6);
	end
	for i=N+1 % second derivative
		t_i = t(1:end-1, N+7)';
		N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, N+2, L);
		J(:,N+7) = N_mtki(:,end-2 : end) * b(1:3,N+7);
	end

	figure
	dx = linspace(x_o(alpha+1), x_o(beta+2), L);
	plot(dx, J(:, 1), 'r', 'LineWidth', 1); %J^{3}_{0}
	hold on;
	plot(dx, J(:, 2), 'g', 'LineWidth', 1); %J^{2}_{0}
	hold on;
	plot(dx, J(:, 3), 'c', 'LineWidth', 1); %J^{1}_{0}**
	hold on;
	plot(dx, J(:, 4), 'k', 'LineWidth', 1); %J^{0}_{0}**
	hold on;
	plot(dx, J(:, 5), 'b', 'LineWidth', 1); %J_{1}**
	hold on;
	plot(dx, J(:, 6), 'm', 'LineWidth', 1); %J_{2}
	hold on;
	plot(dx, J(:, 7), 'y', 'LineWidth', 1); %J_{3}
	hold on;
	plot(dx, J(:, 8), 'r', 'LineWidth', 1); %J_{4}]
	hold on;
	plot(dx, J(:, end-4), 'g', 'LineWidth', 1); %J_{N-1}
	hold on;
	plot(dx, J(:, end-3), 'c', 'LineWidth', 1); %J_{N}
	hold on; 
	plot(dx, J(:, end-2), 'k', 'LineWidth', 1); %J^{0}_{N+1}
	hold on; 
	plot(dx, J(:, end-1), 'b', 'LineWidth', 1); %J^{1}_{N+1}
	hold on;
	plot(dx, J(:, end), 'm', 'LineWidth', 1); %J^{2}_{N+1}
	title('Local interpolating molecules');
	grid minor;
	legend('L_{0}^{3}', 'L_{0}^{2}', 'L_{0}^{1}', 'L_{0}', 'L_{1}', 'L_{2}', 'L_{3}', 'L_{4}', 'L_{N-1}', 'L_{N}', 'L_{N+1}', 'L_{N+1}^{1}', 'L_{N+1}^{2}');

end

end