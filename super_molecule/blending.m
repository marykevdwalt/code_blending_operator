function [B] = blending(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create an array B that contains the blending operator B(f,x) = sum_i f(y_i) M*_{i}(x)
% each row corresponds to an evaluation of B at a different x-value between x_o(alpha) and x_o(beta+1)

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% f_o = original function values at sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

B = zeros(L,1);

if strcmp(knots, 'equal')
    [Mstar] = super_molecule(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots);
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
	B(:) = Mstar(:,1:m-1) * der(1:m-1)' + Mstar(:,0+m:N+1+m) * f_o(1:end)' + Mstar(:,end-m+3:end) * der(m:2*m-3)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), 'k*', 'markersize', 5);
	hold on;
%     plot(x_o, zeros(length(x_o)), 'ko', 'markersize', 5, 'markerfacecolor', 'k');
	dx = linspace(x_o(1), x_o(end), L);
	plot(dx, B, 'r', 'linewidth', 2);
    title('Final result');
    %grid minor;
end

if strcmp(knots, 'half')
    [Mstar] = super_molecule(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots);
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
	B(:) = Mstar(:,1:m-1) * der(1:m-1)' + Mstar(:,0+m:N+m) * f_o(1:end)' + Mstar(:,end-m+2:end) * der(m:2*m-2)';
	
	figure
	s1 = y_o(1:end);
	plot(s1, f_o(1:end), 'k*', 'markersize', 5);
	hold on;
	dx = linspace(x_o(1), x_o(end), L);
%     plot(x_o, zeros(length(x_o)), 'ko', 'markersize', 5, 'markerfacecolor', 'k'); hold on;
	plot(dx, B, 'r', 'linewidth', 2);
    title('Final result');
    %grid minor;
end

end