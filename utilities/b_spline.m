function [N_mxk] = b_spline(x, m, N, k, t)

% evaluate the m-th order b-spline N_{m,x,k}(t) (starting at knot x_k) at t, for k = -m+1,... N
% this is defined (eg) in (6.6.11) in Wavelets: A Mathematical Tool for Signal Analysis by C.K.Chui:
% N_{m,x,k}(t) = (x_{k+m} - x_{k}) [x_{k},...,x_{k+m}](.-t)^{m-1}_{+},  k = -m+1,...,N

if m==1
	if k <= -m
		N_mxk = 0;
	elseif k >= N+1
		N_mxk = 0;
	else
		if t < x(k +m)
			N_mxk = 0;
		elseif t > x(k+1 +m)
			N_mxk = 0;
		else
			N_mxk = 1;
		end
	end

else
	if k <= -m
		N_mxk = 0;
	elseif k >= N+1
		N_mxk = 0;
	else
		TDD = divdiff(x(k+m:k+2*m), max(0, x(k+m:k+2*m) - t).^(m-1));
		N_mxk = (x(k + 2*m) - x(k + m)) * TDD(1, end) ;
	end

end

end



