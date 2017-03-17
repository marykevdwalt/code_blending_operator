function [t] = fine_knot_seq_odd(x, x_o, x_ext, m, N)

% x_o is the original knot sequence used to compute the quasi-interpolating molecule M_i

% create t-knot-sequence to compute local interpolatory molecule
% t-knot-sequence consists of old x-knot-sequence plus (m-1) additional knots taken from x

% x = original sampling points with half-points added in between (necessary to create finer knot sequence t)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

%===================================================================================================================

% m odd:
% insert (m-1) additional knots for each i=0,...,N to facilitate computation of local interpolation operator

% t contains the new extended t knot-sequence for each i
% rows: each t knot-sequence consists of original extended x knot-sequence ((N+2m) points) together with (m-1) additional knots
% columns: i=0,...,N
t = zeros(N+3*m-1, N+1);

%if m=3:
%insert (m-1)=2 additional knots in [x_{i-1}, x_{i}], [x_{i+1}, x_{i+2]
%exceptions near left boundary:
%if i=0, insert 2 knots in [x_0, x_1], [x_1, x_2]
%for i=1, insert 1 stacked knot in [x_0, x_1], 1 knot in [x_1, x_2] (reason: we need 3 stacked knots at left boundary of support of J_1)
%similar exceptions at right boundary

for i=0
	for j=-(m-1):0
		t(j+m,1) = x_ext(j+m);
	end
	for j=1
		t(j+m,1) = x(2);
	end
	for j=2
		t(j+m,1) = x_ext(j+m-1);
	end
	for j=3
		t(j+m,1) = x(4);
	end
	for j=4:N+m+(m-1)
		t(j+m,1) = x_ext(j+m-2);
	end
end

for i=1
	for j=-(m-1):0
		t(j+m,2) = x_ext(j+m);
	end
	for j=1
		t(j+m,2) = x_ext(0+m) +1e-8;
	end
	for j=2:3
		t(j+m,2) = x_ext(j+m-1);
	end
	for j=4
		t(j+m,2) = x(6);
	end
	for j=5:N+m+(m-1)
		t(j+m,2) = x_ext(j+m-2);
	end
end

for i=m-1:N+1-m
	for j=-(m-1):i-1
		t(j+m,i+1) = x_ext(j+m);
	end
	for j=i
		t(j+m, i+1) = x(2*i-1+1);
	end
	for j=i+1:i+2
		t(j+m, i+1) = x_ext(j+m-1);
	end
	for j=i+3
		t(j+m, i+1) = x(2*i+3+1);
	end
	for j=i+4:N+m+(m-1)
		t(j+m,i+1) = x_ext(j+m-2);
	end
end

for i=N-1
	for j=-(m-1):N-2
		t(j+m,i+1) = x_ext(j+m);
	end
	for j=N-1
		t(j+m,i+1) = x(end-5);
	end
	for j=N:N+1
		t(j+m,i+1) = x_ext(j+m-1);
	end
	for j=N+2
		t(j+m,i+1) = x_ext(N+1+m)+1e-8;
	end
	for j=N+3:N+m+(m-1)
		t(j+m,i+1) = x_ext(j+m-2);
	end
end

for i=N
	for j=-(m-1):N-1
		t(j+m,i+1) = x_ext(j+m);
	end
	for j=N
		t(j+m,i+1) = x(end-3);
	end
	for j=N+1
		t(j+m,i+1) = x_ext(j+m-1);
	end
	for j=N+2
		t(j+m,i+1) = x(end-1);
	end
	for j=N+3:N+m+(m-1)
		t(j+m,i+1) = x_ext(j+m-2);
	end
end
	

%disp(t)

end
