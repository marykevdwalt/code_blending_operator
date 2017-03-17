function [J] = interpol_molecule_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create a matrix J that contains the molecules J_{i}(x) = sum_k b_{i,k} N_{m,x,i+k-(m-1)}(x)
% each column corresponds to i
% each row is an evaluation of J_{i}(x) at a different x-value between x_o{alpha} and x_o{beta+1}

% y_o = sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

% implement method in Chui & De Villiers: Applications of optimally local
% interpolation to interpolatory approximants and compactly supported
% wavelets
% Mathematics of Computation of the American Mathematical Society (1966)


if strcmp(knots, 'equal')

    [t_i, t_l, t_r] = fine_knot_seq_new(y_o, x_o, x_ext, m, N, knots);
    t_i = [t_i(1) + [-m+1:-1].*10.^(-3), t_i, t_i(end) + [1:m-1].*10.^(-1)]; % add stacked knots
    [b_i, b_l, b_r] = interpol_molecule_coeff_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots);
    i_m = max(1, floor((m-2)/2));
    J = zeros(L, N+2*m-1);
    
    %% interior molecules
    for i = i_m+1 : N+i_m-m+2
        N_mtki = b_spline_evaluate(x_o, t_i, alpha, beta, m, 2*N+1, L);
        J(:,i+m) = N_mtki(:, 2*i-2*i_m +m : 2*i-2*i_m+m-4 +m) * b_i(:,i);
    end
    
    %% LHS boundary molecules
    
    for i= 1 : i_m
        t_l{i+m} = [x_o(1) + [-m+1:-1].*10.^(-3), x_o(1), t_l{i+m}, x_o(end) + [1:m-1].*10.^(-1)];
        N_mtki = b_spline_evaluate(x_o, t_l{i+m}, alpha, beta, m, 2*N+1, L);
        J(:,i+m) = N_mtki(:, 1 +m : m-3 +m) * b_i(:,i);
    end
    
    % take care of boundaries
    NN = zeros(L,m);
    for l=0:m-1        
        N_mtki = b_spline_evaluate(x_o, t_l{l+1}, alpha, beta, m, m, L);
        NN(:,l+1) = N_mtki(:,l+1);
    end
    for l=0:m-1
        J(:,-l+m) = NN * b_l(:,-l+m);
    end
    
    %% RHS boundary molecules

    for i= N+i_m-m+3 : N
        t_r{i-(N+i_m-m+2)} = [x_o(1) + [-m+1:-1].*10.^(-3), t_r{i-(N+i_m-m+2)}, x_o(end), x_o(end) + [1:m-1].*10.^(-1)];
        N_mtki = b_spline_evaluate(x_o, t_r{i-(N+i_m-m+2)}, alpha, beta, m, length(t_r{i-(N+i_m-m+2)})-2*m, L);
        J(:,i+m) = N_mtki(:, end-2*m+4 : end-m) * b_i(:,i);
    end

    % take care of boundaries
    NN = zeros(L,m-1);
    for r=0:m-2
        N_mtki = b_spline_evaluate(x_o, t_r{m-i_m-1+r}, alpha, beta, m, length(t_r{m-i_m-1+r})-2*m, L);
        NN(:,r+1) = N_mtki(:,end-m+2+r);
    end
    for r=0:m-2
        J(:,N+1+r+m) = NN * b_r(:,r+1);
    end
    
% 	dx = linspace(x_o(alpha+1), x_o(beta+2), L);
% 	figure
% 	plot(dx, J(:,1), 'r'); hold on;
% 	plot(dx, J(:,2), 'g'); hold on;
% 	plot(dx, J(:,3), 'b'); hold on;
% 	plot(dx, J(:,4), 'm'); hold on;
% 	plot(dx, J(:,5), 'k'); hold on;
% 	plot(dx, J(:,6), 'c'); hold on;
% 	plot(dx, J(:,7), 'y'); hold on;
% 	plot(dx, J(:,8), 'r'); hold on;
% 	plot(dx, J(:,9), 'g'); hold on;
% 	plot(dx, J(:,10), 'b'); hold on;
% 	plot(dx, J(:,11), 'm'); hold on;
% 	plot(dx, J(:,12), 'k'); hold on;
% 	plot(dx, J(:,13), 'c'); hold on;
% 	plot(dx, J(:,14), 'y'); hold on;
%   plot(dx, J(:,15), 'r'); hold on;
%   plot(dx, J(:,16), 'b'); hold on;
%   plot(dx, J(:,17), 'g'); hold on;
%   plot(dx, J(:,18), 'c'); hold on;
%   plot(dx, J(:,19), 'k'); hold on;
%   plot(dx, J(:,20), 'm'); hold on;
%   plot(dx, J(:,21), 'r'); hold on;
%   plot(dx, J(:,22), 'b'); hold on;
%   plot(dx, J(:,23), 'g'); hold on;
%   plot(dx, J(:,end-6), 'y'); hold on;
%   plot(dx, J(:,end-5), 'r'); hold on;
% 	plot(dx, J(:,end-4), 'y'); hold on;
% 	plot(dx, J(:,end-3), 'b'); hold on;
% 	plot(dx, J(:,end-2), 'r'); hold on;
% 	plot(dx, J(:,end-1), 'k'); hold on;
% 	plot(dx, J(:,end), 'c'); hold on;
% 	plot(x_o, ones(size(x_o)), '*', 'markersize', 3); hold on;
% 	title('Local interpolating molecules');
%  	grid minor;
end

if strcmp(knots, 'half')

    [t_i, t_l, t_r] = fine_knot_seq_new(y_o, x_o, x_ext, m, N, knots);
    t_i = [t_i(1) + [-m+1:-1].*10.^(-3), t_i, t_i(end) + [1:m-1].*10.^(-1)]; % add stacked knots
    [b_i, b_l, b_r] = interpol_molecule_coeff_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots);
    i_m = max(1, floor((m-2)/2));
    J = zeros(L, N+2*m-1);
    
    %% interior molecules
    for i = i_m+1 : N+i_m-m+1
        N_mtki = b_spline_evaluate(y_o, t_i, alpha, beta, m, 2*N-1, L);
        J(:,i+m) = N_mtki(:, 2*i-2*i_m +m : 2*i-2*i_m+m-4 +m) * b_i(:,i);
    end
    
    %% LHS boundary molecules
    
    for i= 1 : i_m
        N_mtki = b_spline_evaluate(y_o, t_l{i+m}, alpha, beta, m, length(t_l{i+m})-2*m, L);
        J(:,i+m) = N_mtki(:, i-1 +m) / function_evaluate(y_o, N_mtki(:, i-1 +m), i+1);
    end
    
    % take care of boundaries
    NN = zeros(L,m);
    for l=0:m-1        
        N_mtki = b_spline_evaluate(y_o, t_l{l+1}, alpha, beta, m, m, L);
        NN(:,l+1) = N_mtki(:,l+1);
    end
    for l=0:m-1
        J(:,-l+m) = NN * b_l(:,-l+m);
    end
    
    %% RHS boundary molecules

    for i= N+i_m-m+2 : N-1
        N_mtki = b_spline_evaluate(y_o, t_r{i-(N+i_m-m+1)}, alpha, beta, m, length(t_r{i-(N+i_m-m+1)})-2*m, L);
        J(:,i+m) = N_mtki(:, i-1 +m) / function_evaluate(y_o, N_mtki(:, i-1 +m), i+1);
    end

    % take care of boundaries
    NN = zeros(L,m);
    for r=0:m-1
        N_mtki = b_spline_evaluate(y_o, t_r{m-i_m-1+r}, alpha, beta, m, length(t_r{m-i_m-1+r})-2*m, L);
        NN(:,r+1) = N_mtki(:,end-m+1+r);
    end
    for r=0:m-1
        J(:,N+r+m) = NN * b_r(:,r+1);
    end
    
% 	dx = linspace(y_o(alpha+1), y_o(beta+2), L);
% 	figure
% 	plot(dx, J(:,1), 'r'); hold on;
% 	plot(dx, J(:,2), 'g'); hold on;
% 	plot(dx, J(:,3), 'b'); hold on;
% 	plot(dx, J(:,4), 'm'); hold on;
% 	plot(dx, J(:,5), 'k'); hold on;
% 	plot(dx, J(:,6), 'c'); hold on;
% 	plot(dx, J(:,7), 'y'); hold on;
% 	plot(dx, J(:,8), 'r'); hold on;
% 	plot(dx, J(:,9), 'g'); hold on;
% 	plot(dx, J(:,10), 'b'); hold on;
% 	plot(dx, J(:,11), 'm'); hold on;
% 	plot(dx, J(:,12), 'k'); hold on;
% 	plot(dx, J(:,13), 'c'); hold on;
% 	plot(dx, J(:,14), 'y'); hold on;
%   plot(dx, J(:,15), 'r'); hold on;
%   plot(dx, J(:,16), 'b'); hold on;
%   plot(dx, J(:,17), 'g'); hold on;
%   plot(dx, J(:,18), 'c'); hold on;
%   plot(dx, J(:,19), 'k'); hold on;
%   plot(dx, J(:,20), 'm'); hold on;
%   plot(dx, J(:,21), 'r'); hold on;
%   plot(dx, J(:,22), 'b'); hold on;
%   plot(dx, J(:,23), 'g'); hold on;
%   plot(dx, J(:,24), 'm'); hold on;
%   plot(dx, J(:,25), 'k'); hold on;
%   plot(dx, J(:,end-6), 'y'); hold on;
%   plot(dx, J(:,end-5), 'r'); hold on;
% 	plot(dx, J(:,end-4), 'y'); hold on;
% 	plot(dx, J(:,end-3), 'b'); hold on;
% 	plot(dx, J(:,end-2), 'r'); hold on;
% 	plot(dx, J(:,end-1), 'k'); hold on;
% 	plot(dx, J(:,end), 'c'); hold on;
% 	plot(x_o, zeros(size(x_o)), 'ko', 'markersize', 5, 'markerfacecolor', 'k'); hold on;
%   plot(y_o, zeros(size(y_o)), 'ko', 'markersize', 5); hold on;
% 	title('Local interpolating molecules');
%  	grid minor;
end

end