function [] = main_blending()

%% SETUP PARAMETERS

% m is order of B-spline
m = 4;

% knots 'equal' means knots = sampling points
knots = 'equal';
% knots 'half' means knots = halfway between sampling points
% knots = 'half';

% data:
[y_o, y_ext, f_o, x_o, x_ext, N] = data(300, m, knots);

if strcmp(knots, 'equal')
    alpha = 0; % interpolation is done on interval [x_o{alpha}, x_o{beta+1}]
    beta = N;
end

if strcmp(knots, 'half')
    alpha = 0; % interpolation is done on interval [x_o{alpha}, x_o{beta+1}]
    beta = N-1;
end

ll = 15;
L = (N+1)*ll+1;

%% MAIN FUNCTION
B = blending(y_o, y_ext, f_o, x_o, x_ext, alpha, beta, m, N, L, knots);
