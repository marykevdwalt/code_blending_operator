function [N_mxk] = b_spline_2(x, m, k, t)
    %% evaluate N_m(x) = M_m(x) using (6.6.12a)
	%% k is the absolute index shown in (6.6.12a)
	%% set the vector x to be x(1)=x_0=0, x(2)=x_1>0... where x_i is the notation in (6.6.12a)

if m == 1
	N_mxk = (t < x(k+m+2)) && (t >= x(k+m+1)) ; 
else
	a = b_spline_2(x, m-1, k, t) ;	
	b = b_spline_2(x, m-1, k+1, t) ;	
	N_mxk = (t - x(k+m+1)) * a / (x(k+m+m) - x(k+m+1)) + (x(k+m+m+1) - t) * b / (x(k+m+m+1) - x(k+m+2)) ;
end

end
