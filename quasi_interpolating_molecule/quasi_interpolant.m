function [Q] = quasi_interpolant(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create an array Q that contains the quasi-interpolant Q(f,x) = sum_i f(y_i) M_i(x)
% each row corresponds to an evaluation of Q at a different x-value between x_o(alpha) and x_o(beta+1)

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% f_o = original function values at sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

Q = zeros(L,1);

[M] = molecule(y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots);

if strcmp(knots, 'equal')
    TDD = cell(1,2*m-3);
    der = zeros(1,2*m-3);
    for ll=1:m-1
        TDD{m-ll} = divdiff(x_o(1:ll+1), f_o(1:ll+1));
        der(m-ll) = TDD{m-ll}(1,end);
    end
    for r=1:m-2
        TDD{m-1+r} = divdiff(x_o(end-r:end), f_o(end-r:end));
        der(m-1+r) = TDD{m-1+r}(1,end);
    end
	Q(:) = M(:,1:m-1) * der(1:m-1)' + M(:,0+m:N+1+m) * f_o(1:end)' + M(:,end-m+3:end) * der(m:2*m-3)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), '*', 'markersize', 5);
	hold on;
	dx = linspace(x_o(1), x_o(end), L);
	plot(dx, Q, 'r', 'linewidth', 2);
end

if strcmp(knots, 'half')
    TDD = cell(1,2*m-2);
    der = zeros(1,2*m-2);
    for ll=1:m-1
        TDD{m-ll} = divdiff(y_o(1:ll+1), f_o(1:ll+1));
        der(m-ll) = TDD{m-ll}(1,end);
    end
    for r=1:m-1
        TDD{m-1+r} = divdiff(y_o(end-r:end), f_o(end-r:end));
        der(m-1+r) = TDD{m-1+r}(1,end);
    end
	Q(:) = M(:,1:m-1) * der(1:m-1)' + M(:,0+m:N+m) * f_o(1:end)' + M(:,end-m+2:end) * der(m:2*m-2)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), '*', 'markersize', 5);
	hold on;
	dx = linspace(x_o(1), x_o(end), L);
	plot(dx, Q, 'r', 'linewidth', 2);
end

end