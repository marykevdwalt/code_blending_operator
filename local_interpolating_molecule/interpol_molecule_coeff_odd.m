function [b] = interpol_molecule_coeff_odd(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L)

% m odd

% we calculate the coefficients in the interpolating molecule's definition: J_{i}(x) = sum_{k=0}^{2m-2} b_{i,k} N_{m,t,i+k-(m-1)}(x), i=0,...,N
% we do this by requiring that J_{i}(y_j) = delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1),...,i+m-1 (that is, inside support of quasi-interpolating molecule M_{i})

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

% insert (m-1) additional knots for each i to facilitate computation of local interpolation molecule
t = fine_knot_seq_odd(x, x_o, x_ext, m, N);

% b is a matrix containing the coefficients b_{i,k} needed in the definition of the molecule J
% each column corresponds to i=0,...,N
% each row corresponds to k=0,...,2m-2
b = zeros(2*m-1,N+1);

% delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1),...,i+(m-1)
delta = [zeros(m-1,1); 1; zeros(m-1,1)];


%for i=0:N
for i=alpha-m+1:beta+m-1

	%printf("i: %d \n", i)
	
	% support of J_{i} = x_i := [x_{i-(m-1)}, ..., x_{i+m}] = [x_i(1), ..., x_i(2*m)]
	x_i = [x_ext(i-(m-1)+m : i+m+m)]; % add m to make indices positive
	%printf("first_index: %d \n", i-(m-1)+m)
	%printf("last_index: %d \n", i+m+m)
	%printf("length(x_i): %d \n", length(x_i))
	
	% interpolation points inside support of J_{i} = y_i := [y_{i-(m-1)}, ..., y_{i+(m-1)}] = [y_i(1), ..., y_i(2*m-1)]
	y_i = [y_ext(i-(m-1)+m : i+(m-1)+m)]; % add m to make indices positive
	
	% t-knot-sequence to compute splines in J_{i} = t_i := [t_{i-(m-1)}, ..., t_{i+m+m-1}]
	t_i = t(i-(m-1)+m : i+2*m-1+m , i+1)';
	%disp(size(t_i))
	
	% create array to use in function_evaluate.m to get spline values at y_i
	% xy_i := [x_{i-(m-1)}, y_{i-(m-1)}, ..., y_{i+(m-1)}, x_{i+m}] = [xy_i(1), ..., xy_i(2*m+1)]
	xy_i = [x_i(1), y_i, x_i(end)];

	
	% import spline matrix for each i with columns N_{m,t,i-(m-1)},...,N_{m,t,i+m-1}
	% for each i, we use a different knot sequence to compute splines (namely t_{i-(m-1)}... t_{i+2*m-1})
	% in each row, the spline is evaluated at a different value in [t_{i-(m-1)},...,t_{i+2m-1}] = x_i = support of J_{i}
	% we do not extend this knot sequence by adding stacked knots; the stacked knots are already incorporated into the t-sequences
	% number of columns = number of splines = m+N = 2m-1, so N=m-1
	
	%N_mtki = b_spline_evaluate(xy_i, t_i, alpha = first index of xy_i, beta = last index of xy_i - 1, m, m-1, L);
	N_mtki = b_spline_evaluate(xy_i, t_i, 0, 2*m-1, m, m-1, L);
	
	%dx = linspace(xy_i(1), xy_i(end), L);
	%plot(dx, N_mtki(:,1), 'g', 'LineWidth', 1);
	%hold on;
	%plot(dx, N_mtki(:,2), 'r', 'LineWidth', 1);
	%hold on;
	%plot(dx, N_mtki(:,3), 'b', 'LineWidth', 1);
	%hold on;
	%plot(dx, N_mtki(:,4), 'k', 'LineWidth', 1);
	%hold on;
	%plot(dx, N_mtki(:,5), 'm', 'LineWidth', 1);
	
	

	% create matrix with splines N_{m,t,i-(m-1)}(x),...,N_{m,t,i+m-1}(x) evaluated at sampling points on molecule's support (y_j, j=i-(m-1),...,i+m-1)
	S = zeros(2*m-1,2*m-1);
	
	for j=2:2*m
		for k=0:2*m-2
			%printf("j: %d \n", j)
			%printf("k: %d \n", k)
			% evaluate splines N_{m,t,i+k-(m-1)}(x), k=0,...,2m-2 at y_j, j=i-(m-1),...,i+(m-1)
			S(j-1,k+1) = function_evaluate(xy_i, N_mtki(:,k+1), j);
		end
	end

	% [S][b_i] = [delta]; therefore, [b_i] = [inv(S)] * [delta]
	
	b(:,i+1) = S \ delta;
	
end

%disp(b)
end