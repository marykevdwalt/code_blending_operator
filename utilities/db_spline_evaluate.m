function [dN_xmk] = db_spline_evaluate(l, x_o, x_ext, alpha, beta, m, N, L)

% evaluate the l-th derivative of the b-spline N_{m,x,k}(t)

% x_ext is the knot sequence to be used (extended knot sequence is generated in db_spline)
% alpha is the index of the first knot in x_o where we will evaluate the b-spline derivatives
% beta+1 is the index of the last knot in x_o where we will evaluate the b-spline derivatives

% dx is the x-values where we will evaluate the b-spline derivatives
%dx = linspace(x_o(alpha+1), x_o(beta+2), L);
dx = linspace(x_o(1), x_o(end), L);
h = dx(2) - dx(1);

%% one way to evaluate first derivative db-spline:
% evaluate db-spline using an already-evaluated N_{m-1,x,k}-table

% x_ext is the knot sequence to be used
% first evaluate N_{m-1,x,k} using b_spline_evaluate.m
% use these to compute N'_{m,x,k}

% N_xm1k = zeros(L, m+N-1);

% for k = -(m-2):1:N
	% for l = 1: L
		% N_xm1k(l, k+m-1) = b_spline(x_ext(2:end-1), m-1, k, dx(l)) ;
	% end
% end

%dN_xmk = zeros(L, m+N-2);

%for k = -(m-2):1:N-1
%	dN_xmk(:,k +m-1) = ((m-1) / (x_ext(k+m-1 +m)-x_ext(k +m))) * N_xm1k(:,k +m-1) - ((m-1) / (x_ext(k+m +m)-x_ext(k+1 +m))) * N_xm1k(:,k+1 +m-1);
%end

%% a "cleaner" way to evaluate l-th derivative db-spline:
% use db_spline which contains the general formula for the l-th derivative

% we create a matrix containing b-spline derivative values
% each column contains a b-spline l-th derivative starting at a different knot (N_{m,x,-m+1+l)}^(l), ..., N_{m,x,N-l}^(l))
% each row evaluates the b-spline l-th derivatives at a different x-value between x_o(alpha+1) and x_o(beta+2)

dN_xmk = zeros(L, m+N);

for k = -(m-1):1:N

	%printf("k: %d \n", k)
	for j = 1: L
		%printf("t: %d \n", dx(j))
		%dN_xmk(j, k +m) = db_spline(l, x_ext, l+1, numel(x_ext)-l, m, N, k, dx(j)) ;
        dN_xmk(j, k +m) = db_spline(l, x_ext, m, N, k, dx(j)) ;
	end

end

% figure
%  plot(dx, dN_xmk(:,1), 'b', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,2), 'r', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,3), 'g', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,4), 'k', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,5), 'y', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,6), 'c', 'LineWidth', 1);
%  hold on;
%  plot(dx, dN_xmk(:,7), 'm', 'LineWidth', 1);
%  hold on;
%  grid minor;
%  legend('N_{-m+1}', 'N_{-m+2}', 'N_{-m+3}', 'N_{-m+4}', 'N_{-m+5}', 'N_{-m+6}', 'N_{-m+7}');
%  title(['Derivative: l = ', num2str(l)]);
% 
% figure
% plot(dx, dN_xmk(:,end-6), 'r', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end-5), 'g', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end-4), 'b', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end-3), 'y', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end-2), 'm', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end-1), 'c', 'LineWidth', 1);
% hold on;
% plot(dx, dN_xmk(:,end), 'k', 'LineWidth', 1);
% hold on;
% grid minor;
% legend('N_{N-6}', 'N_{N-5}', 'N_{N-4}', 'N_{N-3}', 'N_{N-2}', 'N_{N-1}', 'N_N', 'Location','northwest');
% title(['Derivative: l = ', num2str(l)]);

%disp(dN_xmk)


end
