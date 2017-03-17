function [b] = interpol_molecule_coeff_even(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives)

% m even

% we calculate the coefficients in the interpolating molecule's definition: J_{i}(x) = sum_{k=0}^{2m-3} b_{i,k} N_{m,t,i+k-(m-1)}(x), i=0,...,N+1 plus derivatives
% we do this by requiring that J_{i}(x_j} = delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1)+1,...,i+m-1 (that is, inside support of quasi-interpolating molecule M_{i})

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)


% insert (m-2) additional knots for each i to facilitate computation of local interpolation molecule
t = fine_knot_seq_even(x, x_o, x_ext, m, N, boundaries, derivatives);

% b is a matrix containing the coefficients b_{i,k} needed in the definition of the molecule J
% each column corresponds to i=0,...,N+1 plus derivatives
% each row corresponds to k=0,...,2m-3 (plus extra row for longer b's near boundaries)
b = zeros(2*m-1,N+7);


if (strcmp(boundaries, 'left')) && (strcmp(derivatives, 'one'))
	% interior molecules
	for i=m-1:beta+m-1
		delta = [zeros(m-2,1); 1; zeros(m-1,1)]; % delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1)+1,...,i+m-1
		x_i = [x_ext(i-(m-1)+m : i+m+m)]; % support of J_{i} = x_i := [x_{i-(m-1)}, ..., x_{i+m}] = [x_i(1), ..., x_i(2*m)]
		t_i = t(i-(m-1)+m : i+2*m-2+m , i+4)'; % t-knot-sequence to compute splines in J_{i} = t_i := [t_{i-(m-1)}, ..., t_{i+m+m-2}]
		
		% import spline matrix for each i with columns N_{m,t,i-(m-1)},...,N_{m,t,i+m-2}
		% for each i, we use a different knot sequence to compute splines (namely t_{i-(m-1)}... t_{i+2*m-2})
		% in each row, the spline is evaluated at a different value in [t_{i-(m-1)},...,t_{i+2m-2}] = x_i = support of J_{i}
		% we do not extend this knot sequence by adding stacked knots; the stacked knots are already incorporated into the t-sequences
		% number of columns = number of splines = m+N = 2m-2, so N=m-2
		N_mtki = b_spline_evaluate(x_i, t_i, 0, 2*m-2, m, m-2, L); %N_mtki = b_spline_evaluate(x_i, t_i, alpha = first index of x_i, beta = last index of x_i - 1, m, m-2, L);
			
		% create matrix with splines N_{m,t,i-(m-1)}(x),...,N_{m,t,i+m-2}(x) evaluated at interior knots on molecule's support (x_j, j=i-(m-1)+1,...,i+m-1)
		S = zeros(2*m-2,2*m-2);
		for j=2:2*m-1
			for k=0:2*m-3
				S(j-1,k+1) = function_evaluate(x_i, N_mtki(:,k+1), j); % evaluate splines N_{m,t,i+k-(m-1)}(x), k=0,...,2m-3 at x_j, j=i-(m-1)+1,...,i+m-1
			end
		end

		%for j=1:2*m-2
		%	for k=1:2*m-2
		%		S(j,k) = b_spline(t_i, m, N+2, i+k-(m-1), x_o(i-(m-1)+j +1));
		%	end
		%end
		
		b(1:6,i+4) = S \ delta; % [S][b_i] = [delta]; therefore, [b_i] = [inv(S)] * [delta]
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		delta = [1;0;0;0];
		S = zeros(4,4);
		t_i = t(1:end-2 , i+3)';
		for j=1
			for k=1:4
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(2:5,i+3) = S \ delta;
	end

	for i=0
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=1
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m+1, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+1, x_o(j));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=2
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+2, x_o(j+1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one'))
	% interior molecules
	for i=m-1:N-3
		delta = [zeros(m-2,1); 1; zeros(m-1,1)]; % delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1)+1,...,i+m-1
		x_i = [x_ext(i-(m-1)+m : i+m+m)]; % support of J_{i} = x_i := [x_{i-(m-1)}, ..., x_{i+m}] = [x_i(1), ..., x_i(2*m)]
		t_i = t(i-(m-1)+m : i+2*m-2+m , i+4)'; % t-knot-sequence to compute splines in J_{i} = t_i := [t_{i-(m-1)}, ..., t_{i+m+m-2}]
		
		% import spline matrix for each i with columns N_{m,t,i-(m-1)},...,N_{m,t,i+m-2}
		% for each i, we use a different knot sequence to compute splines (namely t_{i-(m-1)}... t_{i+2*m-2})
		% in each row, the spline is evaluated at a different value in [t_{i-(m-1)},...,t_{i+2m-2}] = x_i = support of J_{i}
		% we do not extend this knot sequence by adding stacked knots; the stacked knots are already incorporated into the t-sequences
		% number of columns = number of splines = m+N = 2m-2, so N=m-2
		N_mtki = b_spline_evaluate(x_i, t_i, 0, 2*m-2, m, m-2, L); %N_mtki = b_spline_evaluate(x_i, t_i, alpha = first index of x_i, beta = last index of x_i - 1, m, m-2, L);
			
		% create matrix with splines N_{m,t,i-(m-1)}(x),...,N_{m,t,i+m-2}(x) evaluated at interior knots on molecule's support (x_j, j=i-(m-1)+1,...,i+m-1)
		S = zeros(2*m-2,2*m-2);
		for j=2:2*m-1
			for k=0:2*m-3
				S(j-1,k+1) = function_evaluate(x_i, N_mtki(:,k+1), j); % evaluate splines N_{m,t,i+k-(m-1)}(x), k=0,...,2m-3 at x_j, j=i-(m-1)+1,...,i+m-1
			end
		end

		%for j=1:2*m-2
		%	for k=1:2*m-2
		%		S(j,k) = b_spline(t_i, m, N+2, i+k-(m-1), x_o(i-(m-1)+j +1));
		%	end
		%end
		
		b(1:6,i+4) = S \ delta; % [S][b_i] = [delta]; therefore, [b_i] = [inv(S)] * [delta]
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		delta = [1;0;0;0];
		S = zeros(4,4);
		t_i = t(1:end-2 , i+3)';
		for j=1
			for k=1:4
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(2:5,i+3) = S \ delta;
	end

	for i=0
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=1
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m+1, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+1, x_o(j));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=2
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+2, x_o(j+1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	% boundary molecules (RIGHT)
	for i=N-2
		delta = [0;0;1;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, N-2+k-3 -1, x_o(N-5+j +1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=N-1
		delta = [0;0;1;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:4
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, N-1+k-3 -1, x_o(N-4+j +1));
			end
		end
		for j=5
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, N-1+k-3 -1, x_o(end));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=N
		delta = [0;0;1;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:4
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, N+k-3 -1, x_o(N-3+j +1));
			end
		end
		for j=5
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, N+k-3 -1, x_o(end));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=N+1
		delta = [0;0;1;0];
		S = zeros(4,4);
		t_i = t(1:end-2 , i+4)';
		for j=1:3
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+1, N+1+k-3 -1, x_o(N-2+j +1));
			end
		end
		for j=4
			for k=1:4
				S(j,k) = db_spline(1, t_i, m, N+1, N+1+k-3 -1, x_o(end));
			end
		end
		b(1:4,i+4) = S \ delta;
	end

	for i=N+1 % first derivative
		delta = [0;0;1];
		S = zeros(3,3);
		t_i = t(1:end-2 , i+5)';
		for j=1:2
			for k=1:3
				S(j,k) = b_spline(t_i, m, N+1, N+2+k-3 -1, x_o(N-1+j +1));
			end
		end
		for j=3
			for k=1:3
				S(j,k) = db_spline(1, t_i, m, N+1, N+2+k-3 -1, x_o(end));
			end
		end
		b(1:3,i+4) = S \ delta;
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one_zero'))
	% interior molecules
	for i=m-1:N-3
		delta = [zeros(m-2,1); 1; zeros(m-1,1)]; % delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1)+1,...,i+m-1
		x_i = [x_ext(i-(m-1)+m : i+m+m)]; % support of J_{i} = x_i := [x_{i-(m-1)}, ..., x_{i+m}] = [x_i(1), ..., x_i(2*m)]
		t_i = t(i-(m-1)+m : i+2*m-2+m , i+4)'; % t-knot-sequence to compute splines in J_{i} = t_i := [t_{i-(m-1)}, ..., t_{i+m+m-2}]
		
		% import spline matrix for each i with columns N_{m,t,i-(m-1)},...,N_{m,t,i+m-2}
		% for each i, we use a different knot sequence to compute splines (namely t_{i-(m-1)}... t_{i+2*m-2})
		% in each row, the spline is evaluated at a different value in [t_{i-(m-1)},...,t_{i+2m-2}] = x_i = support of J_{i}
		% we do not extend this knot sequence by adding stacked knots; the stacked knots are already incorporated into the t-sequences
		% number of columns = number of splines = m+N = 2m-2, so N=m-2
		N_mtki = b_spline_evaluate(x_i, t_i, 0, 2*m-2, m, m-2, L); %N_mtki = b_spline_evaluate(x_i, t_i, alpha = first index of x_i, beta = last index of x_i - 1, m, m-2, L);
			
		% create matrix with splines N_{m,t,i-(m-1)}(x),...,N_{m,t,i+m-2}(x) evaluated at interior knots on molecule's support (x_j, j=i-(m-1)+1,...,i+m-1)
		S = zeros(2*m-2,2*m-2);
		for j=2:2*m-1
			for k=0:2*m-3
				S(j-1,k+1) = function_evaluate(x_i, N_mtki(:,k+1), j); % evaluate splines N_{m,t,i+k-(m-1)}(x), k=0,...,2m-3 at x_j, j=i-(m-1)+1,...,i+m-1
			end
		end

		%for j=1:2*m-2
		%	for k=1:2*m-2
		%		S(j,k) = b_spline(t_i, m, N+2, i+k-(m-1), x_o(i-(m-1)+j +1));
		%	end
		%end
		
		b(1:6,i+4) = S \ delta; % [S][b_i] = [delta]; therefore, [b_i] = [inv(S)] * [delta]
	end
	
	% boundary molecules (LEFT)
	for i=0 % first derivative
		delta = [1;0;0;0];
		S = zeros(4,4);
		t_i = t(1:end-2 , i+3)';
		for j=1
			for k=1:4
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(2:5,i+3) = S \ delta;
	end

	for i=0
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m, x_o(j-1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=1
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1
			for k=1:5
				S(j,k) = db_spline(1, t_i, m, N+1, k-m+1, x_o(1));
			end
		end
		for j=2:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+1, x_o(j));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=2
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, k-m+2, x_o(j+1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	% boundary molecules (RIGHT)
	for i=N-2
		delta = [0;0;1;0;0];
		S = zeros(5,5);
		t_i = t(1:end-2 , i+4)';
		for j=1:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+1, N-2+k-3 -1, x_o(N-5+j +1));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=N-1
		delta = [0;0;1;0];
		S = zeros(4,4);
		t_i = t(1:end-3 , i+4)';
		for j=1:4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N, N-1+k-3 -1, x_o(N-4+j +1));
			end
		end
		b(1:4,i+4) = S \ delta;
	end

	for i=N
		delta = [0;0;1;0];
		S = zeros(4,4);
		t_i = t(1:end-3 , i+4)';
		for j=1:4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N, N+k-3 -1, x_o(N-3+j +1));
			end
		end
		b(1:4,i+4) = S \ delta;
	end

	for i=N+1
		delta = [0;0;1];
		S = zeros(3,3);
		t_i = t(1:end-3 , i+4)';
		for j=1:3
			for k=1:3
				S(j,k) = b_spline(t_i, m, N, N+1+k-3 -1, x_o(N-2+j +1));
			end
		end
		b(1:3,i+4) = S \ delta;
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'three_two'))
	% interior molecules
	for i=m:N-3
		delta = [zeros(m-2,1); 1; zeros(m-1,1)]; % delta_{i,j} = {1 if i=j, \\ 0 otherwise}, j=i-(m-1)+1,...,i+m-1
		x_i = [x_ext(i-(m-1)+m : i+m+m)]; % support of J_{i} = x_i := [x_{i-(m-1)}, ..., x_{i+m}] = [x_i(1), ..., x_i(2*m)]
		t_i = t(i-(m-1)+m : i+2*m-2+m , i+4)'; % t-knot-sequence to compute splines in J_{i} = t_i := [t_{i-(m-1)}, ..., t_{i+m+m-2}]
		
		% import spline matrix for each i with columns N_{m,t,i-(m-1)},...,N_{m,t,i+m-2}
		% for each i, we use a different knot sequence to compute splines (namely t_{i-(m-1)}... t_{i+2*m-2})
		% in each row, the spline is evaluated at a different value in [t_{i-(m-1)},...,t_{i+2m-2}] = x_i = support of J_{i}
		% we do not extend this knot sequence by adding stacked knots; the stacked knots are already incorporated into the t-sequences
		% number of columns = number of splines = m+N = 2m-2, so N=m-2
		N_mtki = b_spline_evaluate(x_i, t_i, 0, 2*m-2, m, m-2, L); %N_mtki = b_spline_evaluate(x_i, t_i, alpha = first index of x_i, beta = last index of x_i - 1, m, m-2, L);
			
		% create matrix with splines N_{m,t,i-(m-1)}(x),...,N_{m,t,i+m-2}(x) evaluated at interior knots on molecule's support (x_j, j=i-(m-1)+1,...,i+m-1)
		S = zeros(2*m-2,2*m-2);
		for j=2:2*m-1
			for k=0:2*m-3
				S(j-1,k+1) = function_evaluate(x_i, N_mtki(:,k+1), j); % evaluate splines N_{m,t,i+k-(m-1)}(x), k=0,...,2m-3 at x_j, j=i-(m-1)+1,...,i+m-1
			end
		end

		%for j=1:2*m-2
		%	for k=1:2*m-2
		%		S(j,k) = b_spline(t_i, m, N+2, i+k-(m-1), x_o(i-(m-1)+j +1));
		%	end
		%end
		
		b(1:6,i+4) = S \ delta; % [S][b_i] = [delta]; therefore, [b_i] = [inv(S)] * [delta]
	end

	% boundary molecules (LEFT)
	for i=0 %third derivative molecule
		delta = [1;0;0;0];
		S = zeros(4,4);
		t_i = t(1:end , 1)';
		for j=1:3
			for k=1:4
				S(j,k) = db_spline(4-j, t_i, m, N+3, k-m, x_o(1));
			end
		end
		for j=4
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+3, k-m, x_o(j-3));
			end
		end
		b(1:4,1) = S \ delta;
	end

	for i=0 %second derivative molecule
		delta = [0;1;0;0;0];
		S = zeros(5,5);
		t_i = t(1:end , 2)';
		for j=1:3
			for k=1:5
				S(j,k) = db_spline(4-j, t_i, m, N+3, k-m, x_o(1));
			end
		end
		for j=4:5
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+3, k-m, x_o(j-3));
			end
		end
		b(1:5,2) = S \ delta;
	end

	for i=0 %first derivative molecule
		delta = [0;0;1;0;0;0];
		S = zeros(6,6);
		t_i = t(1:end , 3)';
		for j=1:3
			for k=1:6
				S(j,k) = db_spline(4-j, t_i, m, N+3, k-m, x_o(1));
			end
		end
		for j=4:6
			for k=1:6
				S(j,k) = b_spline(t_i, m, N+3, k-m, x_o(j-3));
			end
		end
		b(1:6,3) = S \ delta;
	end

	for i=0
		delta = [0;0;0;1;0;0;0];
		S = zeros(7,7);
		t_i = t(1:end , i+4)';
		for j=1:3
			for k=1:7
				S(j,k) = db_spline(4-j, t_i, m, N+3, k-m, x_o(1));
			end
		end
		for j=4:7
			for k=1:7
				S(j,k) = b_spline(t_i, m, N+3, k-m, x_o(j-3));
			end
		end
		b(1:7,i+4) = S \ delta;
	end

	for i=1
		delta = [0;0;0;1;0;0;0];
		S = zeros(7,7);
		t_i = t(1:end , i+4)';
		for j=1:3
			for k=1:7
				S(j,k) = db_spline(4-j, t_i, m, N+3, k-m+1, x_o(1));
			end
		end
		for j=4:7
			for k=1:7
				S(j,k) = b_spline(t_i, m, N+3, k-m+1, x_o(j-2));
			end
		end
		b(1:7,i+4) = S \ delta;
	end

	for i=2
		delta = [0;0;0;1;0;0;0];
		S = zeros(7,7);
		t_i = t(1:end , i+4)';
		for j=1
			for k=1:7
				S(j,k) = db_spline(3, t_i, m, N+3, k-m+2, x_o(1));
			end
		end
		for j=2
			for k=1:7
				S(j,k) = db_spline(2, t_i, m, N+3, k-m+2, x_o(1));
			end
		end
		for j=3:7
			for k=1:7
				S(j,k) = b_spline(t_i, m, N+3, k-m+2, x_o(j-1));
			end
		end
		b(1:7,i+4) = S \ delta;
	end

	for i=3
		delta = [0;0;0;1;0;0;0];
		S = zeros(7,7);
		t_i = t(1:end , i+4)';
		for j=1
			for k=1:7
				S(j,k) = db_spline(3, t_i, m, N+3, k-1, x_o(1));
			end
		end
		for j=2:7
			for k=1:7
				S(j,k) = b_spline(t_i, m, N+3, k-1, x_o(j));
			end
		end
		b(1:7,i+4) = S \ delta;
	end

	% boundary molecules (RIGHT)
	for i=N-2
		delta = [0;0;1;0;0;0];
		S = zeros(6,6);
		t_i = t(1:end-1 , i+4)';
		for j=1:5
			for k=1:6
				S(j,k) = b_spline(t_i, m, N+2, N-2+k-3 -1, x_o(N-5+j +1));
			end
		end
		for j=6
			for k=1:6
				S(j,k) = db_spline(2, t_i, m, N+2, N-2+k-3 -1, x_o(end));
			end
		end
		b(1:6,i+4) = S \ delta;
	end

	for i=N-1
		delta = [0;0;1;0;0;0];
		S = zeros(6,6);
		t_i = t(1:end-1 , i+4)';
		for j=1:4
			for k=1:6
				S(j,k) = b_spline(t_i, m, N+2, N-2+k-3, x_o(N-4+j +1));
			end
		end
		for j=5:6
			for k=1:6
				S(j,k) = db_spline(j-4, t_i, m, N+2, N-2+k-3, x_o(end));
			end
		end
		b(1:6,i+4) = S \ delta;
	end

	for i=N
		delta = [0;0;1;0;0;0];
		S = zeros(6,6);
		t_i = t(1:end-1 , i+4)';
		for j=1:4
			for k=1:6
				S(j,k) = b_spline(t_i, m, N+2, N-4+k, x_o(N-3+j +1));
			end
		end
		for j=5:6
			for k=1:6
				S(j,k) = db_spline(j-4, t_i, m, N+2, N-4+k, x_o(end));
			end
		end
		b(1:6,i+4) = S \ delta;
	end

	for i=N+1
		delta = [0;0;1;0;0];
		S = zeros(5,5);
		t_i = t(1:end-1 , i+4)';
		for j=1:3
			for k=1:5
				S(j,k) = b_spline(t_i, m, N+2, N-3+k, x_o(N-2+j +1));
			end
		end
		for j=4:5
			for k=1:5
				S(j,k) = db_spline(j-3, t_i, m, N+2, N-3+k, x_o(end));
			end
		end
		b(1:5,i+4) = S \ delta;
	end

	for i=N+1 % first derivative molecule
		delta = [0;0;1;0];
		S = zeros(4,4);
		t_i = t(1:end-1 , N+6)';
		for j=1:2
			for k=1:4
				S(j,k) = b_spline(t_i, m, N+2, N-2+k, x_o(N-1+j +1));
			end
		end
		for j=3:4
			for k=1:4
				S(j,k) = db_spline(j-2, t_i, m, N+2, N-2+k, x_o(end));
			end
		end
		b(1:4,N+6) = S \ delta;
	end

	for i=N+1 % second derivative molecule
		delta = [0;0;1];
		S = zeros(3,3);
		t_i = t(1:end-1 , N+7)';
		for j=1
			for k=1:3
				S(j,k) = b_spline(t_i, m, N+2, N-1+k, x_o(N+j +1));
			end
		end
		for j=2:3
			for k=1:3
				S(j,k) = db_spline(j-1, t_i, m, N+2, N-1+k, x_o(end));
			end
		end
		b(1:3,N+7) = S \ delta;
	end

end

end