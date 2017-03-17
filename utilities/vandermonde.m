function [V] = vandermonde(v,r)

% create a Vandermonde matrix with r rows given an array v

n = length(v);

V = zeros(r,n);

for i=1:1:r
	for j=1:1:n
		V(i,j) = (v(j))^(i-1);
	end
end

end