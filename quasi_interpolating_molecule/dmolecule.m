function [dM] = dmolecule(l,y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create a matrix dM that contains the l^th derivatives of the molecules M_{i}(x) = sum_{j=0}^{m-1} a_{i,j} N_{m,x,i+j-(m-1)}(x)
% each column corresponds to i
% each row is an evaluation of M_{i} ' (x) at a different x-value between x_o{alpha} and x_o{beta+1}

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

dM = zeros(L, N+2*m-1);

[a] = molecule_coeff(y_ext, x_ext, alpha, beta, m, N, knots);
dN_mxk = db_spline_evaluate(l,x_o, x_ext, alpha, beta, m, N, L);

if strcmp(knots, 'equal')
	% interior molecules
	for i=0:N
		dM(:,i+m) = dN_mxk(:,i-(m-1)+m : i+m) * a(:,i+m);
	end
	% boundary molecules
    for ll=1:m-1
        dM(:,m-ll) = dN_mxk(:,-(m-1)+m : -ll+m) * a(0+1:m-1-ll+1, m-ll);
    end
    for r=0:m-2
        dM(:,N+1+r+m) = dN_mxk(:,N+r-m+2+m : N+m) * a(r+1:m-2+1, N+1+r+m);
    end
end

if strcmp(knots, 'half')
	% interior molecules
	for i=0:N
		dM(:,i+m) = dN_mxk(:,i-(m-1)+m : i+m) * a(:,i+m);
	end
	% boundary molecules
    for ll=1:m-1
        dM(:,m-ll) = dN_mxk(:,-(m-1)+m : -ll+m) * a(0+1:m-1-ll+1, m-ll);
    end
    for r=1:m-1
        dM(:,N+r+m) = dN_mxk(:,N+r-m+1+m : N+m) * a(r+1:m-1+1, N+r+m);
    end
end

end