function [I] = interpolant(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create an array I that contains the interpolant I(f,x) = sum_i f(y_i) J_i(x)
% each row corresponds to an evaluation of Q at a different x-value between x_o(alpha) and x_o(beta+1)

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% f_o = original function values at sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

I = zeros(L,1);

if strcmp(knots, 'equal')
    [J] = interpol_molecule_new(x, x_o, x_ext, alpha, beta, m, N, L, knots);
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
	I(:) = J(:,1:m-1) * der(1:m-1)' + J(:,0+m:N+1+m) * f_o(1:end)' + J(:,end-m+3:end) * der(m:2*m-3)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), '*', 'markersize', 5);
	hold on;
	dx = linspace(x_o(1), x_o(end), L);
	plot(dx, I, 'r', 'linewidth', 2);
end

if strcmp(knots, 'half')
    [J] = interpol_molecule_new(x, x_o, x_ext, alpha, beta, m, N, L, knots);
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
	I(:) = J(:,1:m-1) * der(1:m-1)' + J(:,0+m:N+m) * f_o(1:end)' + J(:,end-m+2:end) * der(m:2*m-2)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), '*', 'markersize', 5);
	hold on;
	dx = linspace(x_o(1), x_o(end), L);
	plot(dx, I, 'r', 'linewidth', 2);
end
 
end