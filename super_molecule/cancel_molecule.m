function [K] = cancel_molecule(J,y_o, y_ext, x_o, x_ext, alpha, beta, m, N, L, knots)

% we create a matrix K that contains the molecules K_{i}(x) = sum_j M_{i}(y_j) J_{j}(x)
% each column corresponds to i
% each row is an evaluation of K_{i}(x) at a different x-value between x_o{alpha} and x_o{beta+1}

% y_o = sampling points
% y_ext = extended sampling points (with m-1 stacked points at boundaries)
% x_o = original knot sequence
% x_ext = extended knot sequence (with m-1 stacked knots at boundaries)

K = zeros(L, N+2*m-1);

if strcmp(knots, 'equal')
	% get y-values where we need to evaluate M_{i}
    y_j = [y_o(1 : end)];
	% evaluate M_{i} on y_j
	M = molecule(y_o, y_ext, y_j, x_ext, 0, length(y_j)-2, m, N, L, knots);
	% import M_{i}'
    dM = cell(1:m-1);
    for i=1:m-1
        dM_i = dmolecule(i,y_o, y_ext, y_j, x_ext, 0, length(y_j)-2, m, N, L, knots);
        dM{i} = dM_i;
    end
    %% iteration for interior molecules
    % evaluate M_{i} at every y-value in y_j
	M_i = zeros(length(y_j), 1);
    for i = 0:N+1
        for j=1:length(y_j)
			M_i(j) = function_evaluate(y_j, M(:, i+m), j);
        end
        sum_ll = 0;
        sum_r = 0;
        for ll=1:m-1
            tmp = dM{ll}(1,i+m) * J(:,m-ll);
            sum_ll = sum_ll + tmp;
        end
        for r=1:m-2
            tmp = dM{r}(end,i+m) * J(:,end-m+2+r);
            sum_r = sum_r + tmp;
        end
        % construct cancelling molecules (interior)
		K(:,i+m) = J(:, m : end-m+2) * M_i + sum_ll + sum_r;
    end
    %% iteration for boundaries
    % evaluate LHS boundary M_{i} at every y-value in y_j
    M_i_lhs = cell(1,m-1);
    for ll=1:m-1
        M_i_lhs{m-ll} = zeros(length(y_j), 1);
        for j=1:length(y_j)
            M_i_lhs{m-ll}(j) = function_evaluate(y_j, M(:, m-ll), j);
        end
    end
	% evaluate RHS boundary M_{i} at every y-value in y_j
    M_i_rhs = cell(1,m-2);
    for r=1:m-2
        M_i_rhs{r} = zeros(length(y_j), 1);
        for j=1:length(y_j)
            M_i_rhs{r}(j) = function_evaluate(y_j, M(:, end-m+2+r), j);
        end
    end
    % construct cancelling molecule (boundaries)
    % LHS
    for i=1:m-1
        sum_ll = 0;
        for ll=1:m-1
            tmp = dM{ll}(1,m-i) * J(:,m-ll);
            sum_ll = sum_ll + tmp;
        end
        % construct cancelling molecules (interior)
        K(:,m-i) = J(:, m : end-m+2) * M_i_lhs{m-i} + sum_ll;
    end
    % RHS
    for i=N+2:N+m-1
        sum_r = 0;
        for r=1:m-2
            tmp = dM{r}(end,i+m) * J(:,end-m+2+r);
            sum_r = sum_r + tmp;
        end
        % construct cancelling molecules (interior)
        K(:,i+m) = J(:, m : end-m+2) * M_i_rhs{i-(N+1)} + sum_r;
    end
    
end

if strcmp(knots, 'half')
	y_j = [y_o(1 : end)];
	% evaluate M_{i} on y_j
	M = molecule(y_o, y_ext, y_j, x_ext, 0, length(y_j)-2, m, N, L, knots);
	% import M_{i}'
    dM = cell(1:m-1);
    for i=1:m-1
        dM_i = dmolecule(i,y_o, y_ext, y_j, x_ext, 0, length(y_j)-2, m, N, L, knots);
        dM{i} = dM_i;
    end
    %% iteration for interior molecules
    % evaluate M_{i} at every y-value in y_j
	M_i = zeros(length(y_j), 1);
    for i = 0:N
        for j=1:length(y_j)
			M_i(j) = function_evaluate(y_j, M(:, i+m), j);
        end
        sum_ll = 0;
        sum_r = 0;
        for ll=1:m-1
            tmp = dM{ll}(1,i+m) * J(:,m-ll);
            sum_ll = sum_ll + tmp;
        end
        for r=1:m-1
            tmp = dM{r}(end,i+m) * J(:,end-m+1+r);
            sum_r = sum_r + tmp;
        end
        % construct cancelling molecules (interior)
		K(:,i+m) = J(:, m : N +m) * M_i + sum_ll + sum_r;
    end
    %% iteration for boundaries
    % evaluate LHS boundary M_{i} at every y-value in y_j
    M_i_lhs = cell(1,m-1);
    for ll=1:m-1
        M_i_lhs{m-ll} = zeros(length(y_j), 1);
        for j=1:length(y_j)
            M_i_lhs{m-ll}(j) = function_evaluate(y_j, M(:, m-ll), j);
        end
    end
	% evaluate RHS boundary M_{i} at every y-value in y_j
    M_i_rhs = cell(1,m-1);
    for r=1:m-1
        M_i_rhs{r} = zeros(length(y_j), 1);
        for j=1:length(y_j)
            M_i_rhs{r}(j) = function_evaluate(y_j, M(:, end-m+1+r), j);
        end
    end
    % construct cancelling molecule (boundaries)
    % LHS
    for i=1:m-1
        sum_ll = 0;
        for ll=1:m-1
            tmp = dM{ll}(1,m-i) * J(:,m-ll);
            sum_ll = sum_ll + tmp;
        end
        % construct cancelling molecules (interior)
        K(:,m-i) = J(:, m : N +m) * M_i_lhs{m-i} + sum_ll;
    end
    % RHS
    for i=N+1:N+m-1
        sum_r = 0;
        for r=1:m-1
            tmp = dM{r}(end,i+m) * J(:,end-m+1+r);
            sum_r = sum_r + tmp;
        end
        % construct cancelling molecules (interior)
        K(:,i+m) = J(:, m : N +m) * M_i_rhs{i-N} + sum_r;
    end
    
end

end