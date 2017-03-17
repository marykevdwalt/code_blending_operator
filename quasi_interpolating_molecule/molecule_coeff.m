function [a] = molecule_coeff(y_ext, x_ext, alpha, beta, m, N, knots)

% a is a matrix containing the coefficients a_{i,j} needed in the definition of the molecule M
% each column corresponds to i
% each row corresponds to j=0,...,m-1

% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

a = zeros(m,N+2 +3 +2);

if strcmp(knots, 'equal')
	% interior molecules
	for i = m-1:N+1-m
		for j=0:m-1
			denom_a = vandermonde(x_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
		end
	end
	
	% boundary molecules (LEFT)
    for ll = 0:m-1
        for j = 0:m-1-ll
            sigma = symm_poly(m, x_ext(j-(m-1)+1+m : j+m));
            denom_a = vandermonde_confl(x_ext(-m+1+j +m : j +m),m,1,m-j);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, ll+1) = xi;
			a(j+1,m-ll) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = 1:m-2
        for j = 0:m-2-i
            sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m));
            denom_a = vandermonde_confl(x_ext(-m+1+j+i +m : i+j +m),m,1,m-j-i);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = 1:m-2
        for j = m-1-i:m-1
			denom_a = vandermonde(x_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end

	% boundary molecules (RIGHT)
    for r = 0:m-2
        for j = r:m-2
            sigma = symm_poly(m, x_ext(N+1+j-(m-1)+1+m : N+1+j+m));
            denom_a = vandermonde_confl(x_ext(N+j-m+2 +m : N+1+j +m),m,m-j,m);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j+r) = xi;
			a(j+1,N+1+r+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = N-m+2:N
        for j = N-i+1:m-1
            sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m));
            denom_a = vandermonde_confl(x_ext(i+j-m+1 +m : N+1+i+j-N-1 +m),m,N-i-j+m+1,m);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = N+2-m:N
        for j = 0:N-i
			denom_a = vandermonde(x_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end

end

if strcmp(knots, 'half')
	% interior molecules
	for i = m-1:N+1-m
		for j=0:m-1
			denom_a = vandermonde(y_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
		end
	end
	
	% boundary molecules (LEFT)
    for ll = 0:m-1
        for j = 0:m-1-ll
            sigma = symm_poly(m, x_ext(j-(m-1)+1+m : j+m));
            denom_a = vandermonde_confl(y_ext(-m+1+j +m : j +m),m,1,m-j);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, ll+1) = xi;
			a(j+1,m-ll) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = 1:m-2
        for j = 0:m-2-i
            sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m));
            denom_a = vandermonde_confl(y_ext(-m+1+j+i +m : i+j +m),m,1,m-j-i);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = 1:m-2
        for j = m-1-i:m-1
			denom_a = vandermonde(y_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end

	% boundary molecules (RIGHT)
    for r = 1:m-1
        for j = r:m-1
            sigma = symm_poly(m, x_ext(N+j-(m-1)+1+m : N+j+m));
            denom_a = vandermonde_confl(y_ext(N+j-m+1 +m : N+j +m),m,m-j,m);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j+r) = xi;
			a(j+1,N+r+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = N-m+2:N
        for j = N-i+1:m-1
            sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m));
            denom_a = vandermonde_confl(y_ext(i+j-m+1 +m : N+i+j-N +m),m,N-i-j+m,m);
            num_a = denom_a;
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end
    for i = N+2-m:N
        for j = 0:N-i
			denom_a = vandermonde(y_ext(i+j-(m-1)+m : i+j+m),m); % add m to make indices positive
			num_a = denom_a;
			sigma = symm_poly(m, x_ext(i+j-(m-1)+1+m : i+j+m)); % add m to make indices positive
			xi = zeros(m,1);
			for l = 1:m
				xi(l) = sigma(l)/nchoosek(m-1,l-1);
			end
			num_a(:, m-j) = xi;
			a(j+1,i+m) = (det(num_a)) / (det(denom_a));
        end
    end

end

end