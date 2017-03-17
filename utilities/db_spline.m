function [dN_mxk] = db_spline(l, x, m, N, k, t)

% evaluate the l-th derivative of the m-th order b-spline N_{m,x,k}^(l)(t) (with N_{m,x,k} starting at knot x_k) at t, for k = -m+1,... N
% l <= m-1
% this is defined (eg) in Th.6.4 in Wavelets: A Mathematical Tool for Signal Analysis by C.K.Chui:
% N_{m,x,k}'(t) = (m-1)/(x_{k+m-1}-x_{k}) N_{m-1,x,k}(t) - (m-1)/(x_{k+m}-x_{k+1}), N_{m-1,x,k+1}(t)  k = -m+1,...,N
% N_{m,x,k}^(l)(t) = (m-1)/(x_{k+m-1}-x_{k}) N_{m-1,x,k}^(l-1)(t) - (m-1)/(x_{k+m}-x_{k+1}), N_{m-1,x,k+1}^(l-1)(t)  k = -m+1,...,N

if l==1
	if k <= -m
		dN_mxk = 0;
	elseif k >= N+1
		dN_mxk = 0;
	else
		dN_mxk = ((m-1) / (x(k+m-1 +m) - x(k +m)))*b_spline(x(2:end-1), m-1, N, k, t) - ((m-1) / (x(k+m +m) - x(k+1 +m)))*b_spline(x(2:end-1), m-1, N, k+1, t);
	end
else
	if k <= -m
		dN_mxk = 0;
	elseif k > N
		dN_mxk = 0;
	else
		dN_mxk = ((m-1) / (x(k+m-1 +m) - x(k +m)))*db_spline(l-1, x(2:end-1), m-1, N, k, t) - ((m-1) / (x(k+m +m) - x(k+1 +m)))*db_spline(l-1, x(2:end-1), m-1, N, k+1, t);
	end

end

end