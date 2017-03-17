function [b_i, b_l, b_r] = interpol_molecule_coeff_new(y_o, x_o, x_ext, alpha, beta, m, N, L, knots)

% y_o = sampling points
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

% implement method in Chui & De Villiers: Applications of optimally local
% interpolation to interpolatory approximants and compactly supported
% wavelets
% Mathematics of Computation of the American Mathematical Society (1966)

if strcmp(knots, 'equal')
    [t_i, t_l, t_r] = fine_knot_seq_new(y_o, x_o, x_ext, m, N, knots);
    i_m = max(1, floor((m-2)/2));

    delta = [zeros(i_m-1,1); 1; zeros(m-i_m-3,1)];
    b_i = zeros(m-3,N); % interior
    b_l = zeros(m,m); % left hand side boundary
    b_r = zeros(m-1,m-1); % right hand side boundary

    %% coefficients for interior molecules

    for i= i_m+1 : N+i_m-m+2
        S = zeros(m-3,m-3);
        for j=1:m-3
            for k=1:m-3
                S(j,k) = b_spline_wo_boundaries(t_i, m, 2*i-2*i_m+k-1, x_o(i-i_m+j+1));
            end
        end
        b_i(:,i) = S \ delta;
    end

    %% coefficients for LHS boundary molecules

    % replace first i_m-i+1 knots (that coincide with x_o) with unstacked
    % knots to avoid x=x_o
    for i= 1 : i_m
        S = zeros(m-3,m-3);
        for j=i_m-i+1:m-3
            for k=1:m-3
                S(j,k) = b_spline_wo_boundaries(t_l{i+m}, m, k-1, x_o(i-i_m+j+1));                
            end
        end
        S_short = S(i_m-i+1:end, :);
        delta_short = delta(i_m-i+1:end);
        b_i(:,i) = S_short \ delta_short;
    end

    % take care of derivatives
    S = zeros(m,m);
    for j=1
        for k=1:m
            S(j,k) = b_spline(t_l{k},m,N,-m+k,x_o(1));
        end
    end
    for j=2:m
        for k=1:m
            S(j,k) = db_spline(j-1, t_l{k}, m, N, -m+k, x_o(1));
        end
    end
    delta_l = eye(m);
    for i=1:m
        b_l(:,i) = S \ delta_l(:,i);
    end
    b_l = fliplr(b_l);

    %% coefficients for RHS boundary molecules
    
    % replace i-N-i_m+m-2 knots (that coincide with x_{N+1}) with unstacked
    % knots to avoid x=x_{N+1}
    for i= N+i_m-m+3 : N
         S = zeros(m-3,m-3);
         for j=1:N+i_m-i
            for k=1:m-3
                S(j,k) = b_spline_wo_boundaries(t_r{i-(N+i_m-m+2)}, m, (N+i-i_m+m)-m-(m-3)-1+k, x_o(i-i_m+j+1));
            end
        end
         S_short = S(1:end-(i-N-i_m+m-3), :);
         delta_short = delta(1:end-(i-N-i_m+m-3));
         b_i(:,i) = S_short \ delta_short;
    end

    % take care of derivatives
    S = zeros(m-1,m-1);
    for j=1
        for k=1:m-1
            S(j,k) = b_spline(t_r{m-i_m-2+k},m,length(t_r{m-i_m-2+k})-2*m,N,x_o(end));
        end
    end
    for j=2:m-1
        for k=1:m-1
            S(j,k) = db_spline(j-1, t_r{m-i_m-2+k}, m, length(t_r{m-i_m-2+k})-2*m, N, x_o(end));
        end
    end
    delta_r = eye(m-1);
    for i=1:m-1
        b_r(:,i) = S \ delta_r(:,i);
    end
end

if strcmp(knots, 'half')
    [t_i, t_l, t_r] = fine_knot_seq_new(y_o, x_o, x_ext, m, N, knots);
    i_m = max(1, floor((m-2)/2));

    delta = [zeros(i_m-1,1); 1; zeros(m-i_m-3,1)];
    b_i = zeros(m-3,N-1); % interior
    b_l = zeros(m,m); % LHS boundary
    b_r = zeros(m,m); % RHS boundary

    %% coefficients for interior molecules

    for i= i_m+1 : N+i_m-m+1
        S = zeros(m-3,m-3);
        for j=1:m-3
            for k=1:m-3
                S(j,k) = b_spline_wo_boundaries(t_i, m, 2*i-2*i_m+k-1, y_o(i-i_m+j+1));
            end
        end
        b_i(:,i) = S \ delta;
    end

    %% coefficients for LHS boundary molecules
    
    % near LHS boundary, we simply use b-splines with extra knots as
    % molecules (will set this up when creating molecules)
    for i= 1 : i_m
        b_i(:,i) = 1;
    end

    % take care of boundaries
    S = zeros(m,m);
    for j=1
        for k=1:m
            S(j,k) = b_spline(t_l{k},m,N-1,-m+k,y_o(1));
        end
    end
    for j=2:m
        for k=1:m
            S(j,k) = db_spline(j-1, t_l{k}, m, N-1, -m+k, y_o(1));
        end
    end
    delta_l = eye(m);
    for i=1:m
        b_l(:,i) = S \ delta_l(:,i);
    end
    b_l = fliplr(b_l);

    %% coefficients for RHS boundary molecules
    
    % near RHS boundary, we simply use b-splines with extra knots as
    % molecules (will set this up when we create the molecules)
    for i= N+i_m-m+2 : N-1
         b_i(:,i) = 1;
    end

    % take care of derivatives
    S = zeros(m,m);
    for j=1
        for k=1:m
            S(j,k) = b_spline(t_r{m-i_m-2+k},m,length(t_r{m-i_m-2+k})-2*m,N-1,y_o(end));
        end
    end
    for j=2:m
        for k=1:m
            S(j,k) = db_spline(j-1, t_r{m-i_m-2+k}, m, length(t_r{m-i_m-2+k})-2*m, N-1, y_o(end));
        end
    end
    delta_r = eye(m);
    for i=1:m
        b_r(:,i) = S \ delta_r(:,i);
    end

end

end