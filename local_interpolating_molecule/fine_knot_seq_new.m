function [t_i, t_l, t_r] = fine_knot_seq_new(y_o, x_o, x_ext, m, N, knots)

% create finer knot sequence to implement local interpolation operator

if strcmp(knots, 'equal')
    i_m = max(1, floor((m-2)/2));

    %% interior knot sequence: add one extra knot in each interval
    t_i = zeros(1, 2*length(x_o)-1);
    t_i(1:2:end) = x_o;
    tmp = .5 * [x_o(1:end-1) + x_o(2:end)];
    t_i(2:2:end) = tmp;

    %% knot sequence on LHS
    t_l = cell(1,m + i_m);
    for k=0:m-1
        t_l{m-k} = [x_ext(-m+1 +m : -1 +m), linspace(x_o(1), x_o(2), m+1-k), x_o(3:end)];
    end
    for i=1:i_m
        tmp = linspace(t_i(1),t_i(2),i_m-i+3);
        t_l{i+m} = [tmp(2:end-1), t_i(2:end)];
    end

    %% knot sequence on RHS
    t_r = cell(1,m-1 + m-i_m-2);
    for i=N+i_m-m+3:N
        tmp = linspace(t_i(end-1),t_i(end),i-N-i_m+m);
         t_r{i-(N+i_m-m+2)} = [t_i(1:end-1), tmp(2:end-1)];
    end
    for k=0:m-2
         t_r{m-i_m-1+k} = [x_ext(1:N-1+m), linspace(x_o(end-1), x_o(end), m-k), x_ext(end-m+2 : end)];
    end
end

if strcmp(knots, 'half')
    i_m = max(1, floor((m-2)/2));

    %% interior knot sequence: add one extra knot in each interval
    t_i = zeros(1, 2*length(y_o)-1);
    t_i(1:2:end) = y_o;
    t_i(2:2:end) = x_o(2:end-1);

    %% knot sequence on LHS
    t_l = cell(1,m + i_m);
    for k=0:m-1
        t_l{m-k} = [x_ext(-m+1 +m : -1 +m), linspace(y_o(1), y_o(2), m+1-k), y_o(3:end)];
    end
    for i=1:i_m
        tmp = linspace(y_o(i), y_o(i+2), m+1);
        t_l{i+m} = [y_o(1) + [-m+1:-1].*10.^(-3), y_o(1:i-1), tmp, y_o(i+3:end), y_o(end) + [1:m-1].*10.^(-1)];
    end
    t_l{m+1}(m) = 0.5*(t_l{m+1}(m) + t_l{m+1}(m+1));

    %% knot sequence on RHS
    t_r = cell(1,m + m-i_m-2);
    for i=N+i_m-m+2:N-1
        tmp = linspace(y_o(i), y_o(i+2), m+1); % t_r{1},...,t_r{m-i_m-2}
        t_r{i-(N+i_m-m+1)} = [y_o(1) + [-m+1:-1].*10.^(-3), y_o(1:i-1), tmp, y_o(i+3:end), y_o(end) + [1:m-1].*10.^(-1)];
    end
    for k=0:m-1
        t_r{m-i_m-1+k} = [x_ext(-m+1 +m : -1 +m), y_o(1:end-2), linspace(y_o(end-1), y_o(end), m+1-k), x_ext(end-m+2 : end)];
    end
    t_r{end-m}(end-m+1) = 0.5*(t_r{end-m}(end-m) + t_r{end-m}(end-m+1));
end

end