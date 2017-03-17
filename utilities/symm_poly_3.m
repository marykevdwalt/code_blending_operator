function [sigma] = symm_poly_3(y)

% y is an array containing n elements
% n = 3 for symm_poly_3

n = length(y);

sigma = zeros(n+1,1);

sigma(1) = 1;

sigma(2) = 0;
for i=1:1:n
	cur_total = 0;
	cur_total = y(i);
	sigma(2) = sigma(2) + cur_total;
end

sigma(3) = 0;
for i=1:1:n-1
	for j=i+1:1:n
		cur_total = 0;
		cur_total = y(i)*y(j);
		sigma(3) = sigma(3) + cur_total;
	end
end

sigma(4) = 0;
for i=1:1:n-2
	for j=i+1:1:n-1
		for k=j+1:1:n
			cur_total = 0;
			cur_total = y(i)*y(j)*y(k);
			sigma(4) = sigma(4) + cur_total;
		end
	end
end

end