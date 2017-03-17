function [bl,tl,br,tr] = simpler_interpol_molecule_coeff(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives)

% m = 4

% we calculate the local interpolating molecule coefficients
% this is done by simply inserting 1 knot between every two knots
% through this, we create a finer set of b-splines on the new knot sequence
% these new splines ~N_j have the property that ~N_j = 1 at x_{j+1}, j\in\Z, ~N_j=0 for all other x_j

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'three_two'))
	% boundary molecules (LEFT)

	% create knot sequences for boundaries:
	t_00 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, linspace(x_o(1), x_o(2), 5), x_o(3:end)];
	t_01 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1), x_o(1) + .25*(x_o(2)-x_o(1)), x_o(1) + .5*(x_o(2)-x_o(1)), x_o(2), x_o(3:end)];
	t_02 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1), x_o(1) + .5*(x_o(2)-x_o(1)), x_o(2), x_o(3:end)];
	t_03 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1), x_o(2), x_o(3:10)];
	tl = {t_00, t_01, t_02, t_03};

	% bl is a matrix containing the coefficients b_{i,k} needed in the definition of the left-hand side boundary molecules J
	% each column corresponds to i=-3,...,0
	% each row corresponds to k=0,...,3
	bl = zeros(4,4);

	S = zeros(4,4);
	for j=1
		for k=1:4
			S(j,k) = b_spline(tl{k},m,N,1-k,x_o(1));
		end
	end
	for j=2:4
		for k=1:4
			S(j,k) = db_spline(j-1, tl{k}, m, N, 1-k, x_o(1));
		end
	end
	
    %disp('LHS boundary:');
    %disp('S'); disp(S);
	delta = eye(4);
	for i=1:4
		bl(:,i) = S \ delta(:,i);
	end
	bl = fliplr(bl); % this gives the coefficients in the order for J^{3}_{0}, J^{2}_{0}, J^{1}_{0}, J_0
	
	% boundary molecules (RIGHT)

	t_Np10 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1:end-2), x_o(end-1), x_o(end-1) + .5*(x_o(end)-x_o(end-1)), x_o(end-1) + .75*(x_o(end)-x_o(end-1)), x_o(end), x_o(end) + 1e-3, x_o(end) + 2*1e-3, x_o(end) + 3*1e-3];
	t_Np11 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1:end-2), x_o(end-1), x_o(end-1) + .5*(x_o(end)-x_o(end-1)), x_o(end), x_o(end) + 1e-3, x_o(end) + 2*1e-3, x_o(end) + 3*1e-3];
	t_Np12 = [x_o(1) - 3*1e-3, x_o(1) - 2*1e-3, x_o(1) - 1e-3, x_o(1:end-2), x_o(end-1), x_o(end), x_o(end) + 1e-3, x_o(end) + 2*1e-3, x_o(end) + 3*1e-3];
	tr = {t_Np10, t_Np11, t_Np12};

	% br is a matrix containing the coefficients b_{i,k} needed in the definition of the left-hand side boundary molecules J
	% each column corresponds to i=N+1,...,N+3
	% each row corresponds to k=0,...,2
	br = zeros(3,3);

	S = zeros(3,3);
	for j=1
		%S(j,1) = b_spline(tr{1},m,N+2,N,x_o(end));
		%S(j,2) = b_spline(tr{2},m,N+1,N,x_o(end));
		%S(j,3) = b_spline(tr{3},m,N,N,x_o(end));
		for k=1:3
			%S(j,k) = b_spline(tr{k},m,N+3-k,N,x_o(end));
			S(j,k) = b_spline(tr{k},m,length(tr{k})-8,N,x_o(end));
		end
	end
	for j=2:3
		for k=1:3
			%S(j,k) = db_spline(j-1, tr{k}, m, N+3-k, N, x_o(end));
			S(j,k) = db_spline(j-1, tr{k}, m, length(tr{k})-8, N, x_o(end));
		end
	end
	
	%disp('S'); disp(S);
	delta = eye(3);
	for i=1:3
		br(:,i) = S \ delta(:,i);
	end

end

end