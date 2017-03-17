function [sigma] = symm_poly(m,y)

p=poly(y);
n = length(y);
sigma = zeros(n+1,1);
for k = 1:n+1
    sigma(k)=(-1)^(k+1) * p(k);
end

end