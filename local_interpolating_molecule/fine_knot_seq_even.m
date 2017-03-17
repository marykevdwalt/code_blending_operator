function [t] = fine_knot_seq_even(x, x_o, x_ext, m, N, boundaries, derivatives)

% x_o is the original knot sequence used to compute the quasi-interpolating molecule M_i

% create t-knot-sequence to compute local interpolatory molecule
% t-knot-sequence consists of old x-knot-sequence plus (m-2) additional knots taken from x

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

%===================================================================================================================

% m even:
% insert (m-2) additional knots for each i=m-1,...,N-3 to facilitate computation of local interpolation operator
% exceptions near boundaries (i=0, i=N+1)

% t contains the new extended t knot-sequence for each i
% rows: each t knot-sequence consists of original extended x knot-sequence ((N+2m) points) together with (m-2) additional knots
% columns: i=0,...,N+1 plus 3 columns on LHS and 2 columns on RHS for boundary molecules
t = zeros(N+3*m-1, N+7); %# original knots = N+2m; at most we add m-1 additional knots

%details:
%insert (m-2)=2 additional knots in [x_{i-1}, x_{i}], [x_{i}, x_{i+1}]
%exceptions near left boundary:
%if i=0, third derivative, insert 3 knots in [x_0, x_1]
%if i=0, second derivative, insert 1 knot in [x_0, x_1], 2 knots in [x_1,x_2]
%if i=0, first derivative, insert 1 knot in [x_0, x_1], [x_1,x_2], [x_2,x_3]
%if i=0, insert 1 knot in [x_0, x_1]. [x_1,x_2], [x_2,x_3]
%if i=1, insert 1 knot in [x_0, x_1], [x_1,x_2], [x_2,x_3]
%if i=2, insert 1 knot in [x_1, x_2], [x_2,x_3], [x_3,x_4]
%if i=3, insert 1 knot in [x_2, x_3], [x_3,x_4], [x_4,x_5]
%exceptions near right boundary:
%if i=N-2, insert 1 knot in [x_N, x_N+1]
%if i=N-1, insert 1 knot in [x_N, x_N+1]
%if i=N, insert 1 knot in [x_N, x_N+1]
%if i=N+1, insert 1 knot in [x_N, x_N+1]
%if i=N+2, insert 1 knot in [x_N, x_N+1]

if (strcmp(boundaries, 'left')) && (strcmp(derivatives, 'one'))
	% boundary molecules (LEFT)
	for i=0 % first derivative
		for j=-(m-1):0
			t(j+m,i+3) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+3) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+3) = x_ext(j+m-1);
		end
	end
	for i=0
		for j=-(m-1):0
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+4) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=1
		for j=-(m-1):0
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+4) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=2
		for j=-(m-1):1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=2
			t(j+m,i+4) = x(4);
		end
		for j=3:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end

	% interior molecules
	for i=m-1:N-3
		for j=-(m-1):i-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=i
			t(j+m, i+4) = x(2*i-1+1);
		end
		for j=i+1
			t(j+m, i+4) = x_ext(j+m-1);
		end
		for j=i+2
			t(j+m, i+4) = x(2*i+1+1);
		end
		for j=i+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

end


if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one'))
	% boundary molecules (LEFT)
	% for i=0 % first derivative % 1 knot in [x_0,x_1] %% ORIGINAL
		% for j=-(m-1):0
			% t(j+m,i+3) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,i+3) = x(2);
		% end
		% for j=2:N+m+(m-2)-1
			% t(j+m,i+3) = x_ext(j+m-1);
		% end
	% end
	for i=0 % first derivative % 1 knot in [x_1,x_2]
		for j=-(m-1):1
			t(j+m,i+3) = x_ext(j+m);
		end
		for j=2
			t(j+m,i+3) = x_o(2) + 0.5*(x_o(3)-x_o(2));
		end
		for j=3:N+m+(m-2)-1
			t(j+m,i+3) = x_ext(j+m-1);
		end
	end
	% for i=0 % 1 knot in [x_0,x_1] %% 	ORIGINAL
		% for j=-(m-1):0
			% t(j+m,i+4) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,i+4) = x(2);
		% end
		% for j=2:N+m+(m-2)-1
			% t(j+m,i+4) = x_ext(j+m-1);
		% end
	% end
	for i=0 % 1 knot in [x_1,x_2]
		for j=-(m-1):1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=2
			t(j+m,i+4) = x_o(2) + 0.5*(x_o(3)-x_o(2));
		end
		for j=3:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=1
		for j=-(m-1):0
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+4) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=2
		for j=-(m-1):1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=2
			t(j+m,i+4) = x(4);
		end
		for j=3:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end

	% interior molecules
	for i=m-1:N-3
		for j=-(m-1):i-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=i
			t(j+m, i+4) = x(2*i-1+1);
		end
		for j=i+1
			t(j+m, i+4) = x_ext(j+m-1);
		end
		for j=i+2
			t(j+m, i+4) = x(2*i+1+1);
		end
		for j=i+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	% boundary molecules (RIGHT)
	for i=N-2
		for j=-(m-1):N
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+4) = x(2*N+1+1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=N-1
		for j=-(m-1):N
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+4) = x(end-1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=N
		for j=-(m-1):N
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+4) = x(end-1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=N+1
		for j=-(m-1):N
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+4) = x(end-1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=N+1 % first derivative
		for j=-(m-1):N
			t(j+m,i+5) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+5) = x(end-1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+5) = x_ext(j+m-1);
		end
	end

end

if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'one_zero'))
	% boundary molecules (LEFT)
	for i=0 % first derivative
		for j=-(m-1):0
			t(j+m,i+3) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+3) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+3) = x_ext(j+m-1);
		end
	end
	for i=0
		for j=-(m-1):0
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+4) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=1
		for j=-(m-1):0
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=1
			t(j+m,i+4) = x(2);
		end
		for j=2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=2
		for j=-(m-1):1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=2
			t(j+m,i+4) = x(4);
		end
		for j=3:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end

	% interior molecules
	for i=m-1:N-3
		for j=-(m-1):i-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=i
			t(j+m, i+4) = x(2*i-1+1);
		end
		for j=i+1
			t(j+m, i+4) = x_ext(j+m-1);
		end
		for j=i+2
			t(j+m, i+4) = x(2*i+1+1);
		end
		for j=i+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	% boundary molecules (RIGHT)
	for i=N-2
		for j=-(m-1):N
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,i+4) = x(2*N+1+1);
		end
		for j=N+2:N+m+(m-2)-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
	end
	for i=N-1
		for j=-(m-1):N+m+(m-2)-2
			t(j+m,i+4) = x_ext(j+m);
		end
	end
	for i=N
		for j=-(m-1):N+m+(m-2)-2
			t(j+m,i+4) = x_ext(j+m);
		end
	end
	for i=N+1
		for j=-(m-1):N+m+(m-2)-2
			t(j+m,i+4) = x_ext(j+m);
		end
	end

end

if (strcmp(boundaries, 'both')) && (strcmp(derivatives, 'three_two'))
	% boundary molecules (LEFT)
	for i=0 %third derivative
		for j=-(m-1):0
			t(j+m,1) = x_ext(j+m);
		end
		for j=1
			t(j+m,1) = x_o(1) + .25*(x_o(2)-x_o(1));
		end
		for j=2
			t(j+m,1) = x_o(1) + .5*(x_o(2)-x_o(1));
		end
		for j=3
			t(j+m,1) = x_o(1) + .75*(x_o(2)-x_o(1));
		end
		for j=4:N+2*m-1
			t(j+m,1) = x_ext(j-3+m);
		end
	end

	for i=0 %second derivative % 1 knot in [x_0,x_1] and 2 knots in [x_1,x_2] %% ORIGINAL
		for j=-(m-1):0
			t(j+m,2) = x_ext(j+m);
		end
		for j=1
			t(j+m,2) = x_o(1)+.5*(x_o(2)-x_o(1));
		end
		for j=2
			t(j+m,2) = x_ext(j+m-1);
		end
		for j=3
			t(j+m,2) = x_o(2) + .33*(x_o(3)-x_o(2));
		end
		for j=4
			t(j+m,2) = x_o(2) + .67*(x_o(3)-x_o(2));
		end
		for j=5:N+2*m-1
			t(j+m,2) = x_ext(j+m-3);
		end
	end
	% for i=0 %second derivative %3 knots in [x_1,x_2]
		% for j=-(m-1):1
			% t(j+m,2) = x_ext(j+m);
		% end
		% for j=2
			% t(j+m,2) = x_o(2)+.25*(x_o(3)-x_o(2));
		% end
		% for j=3
			% t(j+m,2) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=4
			% t(j+m,2) = x_o(2) + .75*(x_o(3)-x_o(2));
		% end
		% for j=5:N+2*m-1
			% t(j+m,2) = x_ext(j+m-3);
		% end
	% end

	for i=0 %first derivative % 1 knot in each of [x_0,x_1], [x_1,x_2], [x_2,x_3] %% ORIGINAL
		for j=-(m-1):0
			t(j+m,3) = x_ext(j+m);
		end
		for j=1
			t(j+m,3) = x_o(1)+.5*(x_o(2)-x_o(1));
		end
		for j=2
			t(j+m,3) = x_ext(j+m-1);
		end
		for j=3
			t(j+m,3) = x_o(2) + .5*(x_o(3)-x_o(2));
		end
		for j=4
			t(j+m,3) = x_ext(j+m-2);
		end
		for j=5
			t(j+m,3) = x_o(3) + .5*(x_o(4)-x_o(3));
		end
		for j=6:N+2*m-1
			t(j+m,3) = x_ext(j+m-3);
		end
	end
	% for i=0 %first derivative % 1 knot in [x_1,x_2], 2 knots in [x_2,x_3] %% BLOWS UP
		% for j=-(m-1):1
			% t(j+m,3) = x_ext(j+m);
		% end
		% for j=2
			% t(j+m,3) = x_o(2)+.5*(x_o(3)-x_o(2));
		% end
		% for j=3
			% t(j+m,3) = x_ext(j+m-1);
		% end
		% for j=4
			% t(j+m,3) = x_o(3)+.33*(x_o(4)-x_o(3));
		% end
		% for j=5
			% t(j+m,3) = x_o(3) + .66*(x_o(4)-x_o(3));
		% end
		% for j=6:N+2*m-1
			% t(j+m,3) = x_ext(j+m-3);
		% end
	% end
	% for i=0 %first derivative % 2 knots in [x_1,x_2], 1 knot in [x_2,x_3] %% BLOWS UP
		% for j=-(m-1):1
			% t(j+m,3) = x_ext(j+m);
		% end
		% for j=2
			% t(j+m,3) = x_o(2)+.33*(x_o(3)-x_o(2));
		% end
		% for j=3
			% t(j+m,3) = x_o(2)+.67*(x_o(3)-x_o(2));
		% end
		% for j=4
			% t(j+m,3) = x_ext(j+m-2);
		% end
		% for j=5
			% t(j+m,3) = x_o(3) + .5*(x_o(4)-x_o(3));
		% end
		% for j=6:N+2*m-1
			% t(j+m,3) = x_ext(j+m-3);
		% end
	% end
	% for i=0 %first derivative % 2 knots in [x_0,x_1], 1 knot in [x_1,x_2]
		% for j=-(m-1):0
			% t(j+m,3) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,3) = x_o(1)+.33*(x_o(2)-x_o(1));
		% end
		% for j=2
			% t(j+m,3) = x_o(1)+.67*(x_o(2)-x_o(1));
		% end
		% for j=3
			% t(j+m,3) = x_ext(j+m-2);
		% end
		% for j=4
			% t(j+m,3) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=5:N+2*m-1
			% t(j+m,3) = x_ext(j+m-3);
		% end
	% end

	% for i=0 % 1 knot in each of [x_0,x_1], [x_1,x_2], [x_2,x_3] %% ORIGINAL
		% for j=-(m-1):0
			% t(j+m,4) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,4) = x_o(1)+.5*(x_o(2)-x_o(1));
		% end
		% for j=2
			% t(j+m,4) = x_ext(j+m-1);
		% end
		% for j=3
			% t(j+m,4) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=4
			% t(j+m,4) = x_ext(j+m-2);
		% end
		% for j=5
			% t(j+m,4) = x_o(3) + .5*(x_o(4)-x_o(3));
		% end
		% for j=6:N+2*m-1
			% t(j+m,4) = x_ext(j+m-3);
		% end
	% end
	for i=0 % 1 knot in each of [x_1,x_2], [x_2,x_3], [x_3,x_4]
		for j=-(m-1):1
			t(j+m,4) = x_ext(j+m);
		end
		for j=2
			t(j+m,4) = x_o(2)+.5*(x_o(3)-x_o(2));
		end
		for j=3
			t(j+m,4) = x_ext(j+m-1);
		end
		for j=4
			t(j+m,4) = x_o(3) + .5*(x_o(4)-x_o(3));
		end
		for j=5
			t(j+m,4) = x_ext(j+m-2);
		end
		for j=6
			t(j+m,4) = x_o(4) + .5*(x_o(5)-x_o(4));
		end
		for j=7:N+2*m-1
			t(j+m,4) = x_ext(j+m-3);
		end
	end
	% for i=0 % 2 knots in [x_1,x_2], 1 knot in [x_2,x_3]
		% for j=-(m-1):1
			% t(j+m,4) = x_ext(j+m);
		% end
		% for j=2
			% t(j+m,4) = x_o(1)+.33*(x_o(2)-x_o(1));
		% end
		% for j=3
			% t(j+m,4) = x_o(1)+.67*(x_o(2)-x_o(1));
		% end
		% for j=4
			% t(j+m,4) = x_ext(j+m-2);
		% end
		% for j=5
			% t(j+m,4) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=6:N+2*m-1
			% t(j+m,4) = x_ext(j+m-3);
		% end
	% end

	% for i=1 % 1 knot in each of [x_0,x_1], [x_1,x_2], [x_2,x_3] %% ORIGINAL
		% for j=-(m-1):0
			% t(j+m,5) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,5) = x_o(1)+.5*(x_o(2)-x_o(1));
		% end
		% for j=2
			% t(j+m,5) = x_ext(j+m-1);
		% end
		% for j=3
			% t(j+m,5) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=4
			% t(j+m,5) = x_ext(j+m-2);
		% end
		% for j=5
			% t(j+m,5) = x_o(3) + .5*(x_o(4)-x_o(3));
		% end
		% for j=6:N+2*m-1
			% t(j+m,5) = x_ext(j+m-3);
		% end
	% end
	for i=1 % 1 knot in each of [x_1,x_2], [x_2,x_3], [x_3,x_4]
		for j=-(m-1):1
			t(j+m,5) = x_ext(j+m);
		end
		for j=2
			t(j+m,5) = x_o(2)+.5*(x_o(3)-x_o(2));
		end
		for j=3
			t(j+m,5) = x_ext(j+m-1);
		end
		for j=4
			t(j+m,5) = x_o(3) + .5*(x_o(4)-x_o(3));
		end
		for j=5
			t(j+m,5) = x_ext(j+m-2);
		end
		for j=6
			t(j+m,5) = x_o(4) + .5*(x_o(5)-x_o(4));
		end
		for j=7:N+2*m-1
			t(j+m,5) = x_ext(j+m-3);
		end
	end
	% for i=1 % 2 knots in [x_0,x_1], 1 knot in [x_1,x_2]
		% for j=-(m-1):0
			% t(j+m,5) = x_ext(j+m);
		% end
		% for j=1
			% t(j+m,5) = x_o(1)+.33*(x_o(2)-x_o(1));
		% end
		% for j=2
			% t(j+m,5) = x_o(1)+.67*(x_o(2)-x_o(1));
		% end
		% for j=3
			% t(j+m,5) = x_ext(j+m-2);
		% end
		% for j=4
			% t(j+m,5) = x_o(2) + .5*(x_o(3)-x_o(2));
		% end
		% for j=5:N+2*m-1
			% t(j+m,5) = x_ext(j+m-3);
		% end
	% end

	for i=2 %% ORIGINAL
		for j=-(m-1):1
			t(j+m,6) = x_ext(j+m);
		end
		for j=2
			t(j+m,6) = x_o(2)+.5*(x_o(3)-x_o(2));
		end
		for j=3
			t(j+m,6) = x_ext(j+m-1);
		end
		for j=4
			t(j+m,6) = x_o(3) + .5*(x_o(4)-x_o(3));
		end
		for j=5
			t(j+m,6) = x_ext(j+m-2);
		end
		for j=6
			t(j+m,6) = x_o(4) + .5*(x_o(5)-x_o(4));
		end
		for j=7:N+2*m-1
			t(j+m,6) = x_ext(j+m-3);
		end
	end

	for i=3 %% ORIGINAL
		for j=-(m-1):2
			t(j+m,7) = x_ext(j+m);
		end
		for j=3
			t(j+m,7) = x_o(3)+.5*(x_o(4)-x_o(3));
		end
		for j=4
			t(j+m,7) = x_ext(j+m-1);
		end
		for j=5
			t(j+m,7) = x_o(4) + .5*(x_o(5)-x_o(4));
		end
		for j=6
			t(j+m,7) = x_ext(j+m-2);
		end
		for j=7
			t(j+m,7) = x_o(5) + .5*(x_o(6)-x_o(5));
		end
		for j=8:N+2*m-1
			t(j+m,7) = x_ext(j+m-3);
		end
	end

	% interior molecules
	for i=m:N-3
		for j=-(m-1):i-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=i
			t(j+m, i+4) = x(2*i-1+1);
		end
		for j=i+1
			t(j+m, i+4) = x_ext(j+m-1);
		end
		for j=i+2
			t(j+m, i+4) = x(2*i+1+1);
		end
		for j=i+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	% boundary molecules (RIGHT)
	for i=N-2 % 1 knot each in [x_{N-3},x_{N-2}], [x_{N-2},x_{N-1}]
		for j=-(m-1):N-3
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N-2
			t(j+m,i+4) = x_o(N-2) + .5*(x_o(N-1)-x_o(N-2));
		end
		for j=N-1
			t(j+m,i+4) = x_ext(j+m-1);
		end
		for j=N
			t(j+m,i+4) = x_o(N-1) + .5*(x_o(N)-x_o(N-1));
		end
		for j=N+1:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	for i=N-1 % 1 knot each in [x_{N-2},x_{N-1}], [x_{N-1},x_{N}]
		for j=-(m-1):N-2
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N-1
			t(j+m,i+4) = x_o(N-1) + .5*(x_o(N)-x_o(N-1));
		end
		for j=N
			t(j+m,i+4) = x_ext(j+m-1);
		end
		for j=N+1
			t(j+m,i+4) = x_o(N) + .5*(x_o(N+1)-x_o(N));
		end
		for j=N+2:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	for i=N % 1 knot each in [x_{N-1},x_{N}], [x_{N},x_{N+1}]
		for j=-(m-1):N-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N
			t(j+m,i+4) = x_o(N) + .5*(x_o(N+1)-x_o(N));
		end
		for j=N+1
			t(j+m,i+4) = x_ext(j+m-1);
		end
		for j=N+2
			t(j+m,i+4) = x_o(N+1) + .5*(x_o(N+2)-x_o(N+1));
		end
		for j=N+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	for i=N+1 % 1 knot each in [x_{N-1},x_{N}], [x_{N},x_{N+1}]
		for j=-(m-1):N-1
			t(j+m,i+4) = x_ext(j+m);
		end
		for j=N
			t(j+m,i+4) = x_o(N) + .5*(x_o(N+1)-x_o(N));
		end
		for j=N+1
			t(j+m,i+4) = x_ext(j+m-1);
		end
		for j=N+2
			t(j+m,i+4) = x_o(N+1) + .5*(x_o(N+2)-x_o(N+1));
		end
		for j=N+3:N+2*m-1-1
			t(j+m,i+4) = x_ext(j+m-2);
		end
	end

	for i=N+1 %first derivative % 1 knot each in [x_{N-1},x_{N}], [x_{N},x_{N+1}]
		for j=-(m-1):N-1
			t(j+m,N+6) = x_ext(j+m);
		end
		for j=N
			t(j+m,N+6) = x_o(N) + .5*(x_o(N+1)-x_o(N));
		end
		for j=N+1
			t(j+m,N+6) = x_ext(j+m-1);
		end
		for j=N+2
			t(j+m,N+6) = x_o(N+1) + .5*(x_o(N+2)-x_o(N+1));
		end
		for j=N+3:N+2*m-1-1
			t(j+m,N+6) = x_ext(j+m-2);
		end
	end

	for i=N+1 %second derivative % 2 knots in [x_{N}, x_{N+1}]
		for j=-(m-1):N
			t(j+m,N+7) = x_ext(j+m);
		end
		for j=N+1
			t(j+m,N+7) = x_o(N+1) + .33*(x_o(N+2)-x_o(N+1));
		end
		for j=N+2
			t(j+m,N+7) = x_ext(j+m-1);
		end
		for j=N+3
			t(j+m,N+7) = x_o(N+1) + .67*(x_o(N+2)-x_o(N+1));
		end
		for j=N+4:N+2*m-1-1
			t(j+m,N+7) = x_ext(j+m-2);
		end
	end

end

end
