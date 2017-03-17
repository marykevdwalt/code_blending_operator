function [CV] = vandermonde_confl(v, r, r_1, r_2)

% create confluent vandermonde matrix with r rows for an array v with v(r_{1}) = v(r_{1}+1) = ... = v(r_{2})
% vandermonde_confl.m calls vandermonde.m

n = length(v);

V = vandermonde(v,r);

CV = V;

% replace repeating columns with k^th derivative (k = 1,...,r_2 - r_1)

for i = 1:1:r
	for k = 1:1:r_2-r_1
		if k >= i
			CV(i, r_1+k) = 0;
		else
			CV(i,r_1+k) = (factorial(i-1) / factorial(i-1-k)) * (v(r_1))^(i-k-1);
		end	
	end
end



end