function [J] = simpler_interpol_molecule(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives)

% we create a matrix J that contains the molecules J_{i}(x)
% each J_{i}(x) is simply a b-spline on the original knot sequence with 1 knot added in between every two knots (the x-sequence, therefore)
% each column corresponds to i=0,...,N+1 plus derivatives
% each row is an evaluation of J_{i}(x) at a different x-value between x_ext(1)=t(1,i+1) and x_ext(end)=t(end,i+1)

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

J = zeros(L, N+7);

if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'three_two'))
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
	
	% interior molecules
	for i=2:N
		%printf('i = %d \n', i);
		%printf('max_value = %d \n', double(max(N_mtki(:,2*(i-1) +m))));
		%printf('function_evaluate = %d \n', double(function_evaluate(x_o, N_mtki(:,2*(i-1) +m), i+1)));
		J(:,i+4) = N_mtki(:,2*(i-1) +m) / function_evaluate(x_o, N_mtki(:,2*(i-1) +m), i+1);
		%J(:,i+4) = N_mtki(:,2*(i-1) +m) / max(N_mtki(:,2*(i-1) +m));
	end
	for i=1
		new_knots = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1), x_o(1) + 0.25*(x_o(2)-x_o(1)), x_o(1) + 0.5*(x_o(2)-x_o(1)), x_o(2), x_o(2) + 0.5*(x_o(3)-x_o(2)), x_o(3), x_o(4:end), x_o(end) + 1e-3, x_o(end) + 2*1e-3, x_o(end) + 3*1e-3];
		NN = b_spline_evaluate(x_o, new_knots, alpha, beta, m, length(new_knots)-8, L); N_1 = NN(:,1 +4);
		J(:,i+4) = N_1 / function_evaluate(x_o, N_1, 2);
	end
	
	% boundary molecules (LEFT)
	for i=-3:0
		J(:,i+4) = N_bdry_left * bl(:,i+4);
	end

	% boundary molecules (RIGHT)
	for i=N+1:N+3
		J(:,i+4) = N_bdry_right * br(:,i-N);
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
	plot(dx, J(:,9), 'g'); hold on;
	plot(dx, J(:,10), 'b'); hold on;
	plot(dx, J(:,11), 'm'); hold on;
	plot(dx, J(:,12), 'k'); hold on;
	plot(dx, J(:,13), 'c'); hold on;
	plot(dx, J(:,14), 'y'); hold on;
    plot(dx, J(:,15), 'r'); hold on;
    plot(dx, J(:,16), 'b'); hold on;
    plot(dx, J(:,17), 'g'); hold on;
    plot(dx, J(:,18), 'c'); hold on;
    plot(dx, J(:,19), 'k'); hold on;
    plot(dx, J(:,20), 'm'); hold on;
    plot(dx, J(:,21), 'r'); hold on;
    plot(dx, J(:,22), 'b'); hold on;
    plot(dx, J(:,23), 'g'); hold on;
    plot(dx, J(:,24), 'k'); hold on;
	%plot(dx, J(:,end-5), 'r'); hold on;
	%plot(dx, J(:,end-4), 'g'); hold on;
	plot(dx, J(:,end-3), 'b'); hold on;
	plot(dx, J(:,end-2), 'm'); hold on;
	plot(dx, J(:,end-1), 'k'); hold on;
	plot(dx, J(:,end), 'c'); hold on;
	plot(x_o, ones(size(x_o)), '*', 'markersize', 3); hold on;
	%legend('J^3_0', 'J^2_0', 'J^1_0', 'J_0', 'J_1', 'J_2', 'J_3', 'J_4', 'J_5', 'J_6', 'J_7', 'J_8', 'J_9', 'J_10', 'J_{N-2}', 'J_{N-1}', 'J_N', 'J_{N+1}', 'J^1_{N+1}', 'J^2_{N+1}');
	title('Local interpolating molecules');
	grid minor;

end

end