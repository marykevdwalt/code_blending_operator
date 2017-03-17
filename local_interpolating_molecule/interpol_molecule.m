function [J] = interpol_molecule(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives)

if m == 3
	J = interpol_molecule_odd(y_o, y_ext, x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives);
	
elseif m == 4
	J = interpol_molecule_even(x, x_o, x_ext, alpha, beta, m, N, L, boundaries, derivatives);

end