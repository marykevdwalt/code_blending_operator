% spline derivatives at boundaries

m=4;
x_o = [0:1:10];
N = length(x_o) - 2;
x_ext = [ x_o(1) + [-m+1:-1].*10.^(-3) , x_o , x_o(end) + [1:m-1].*10.^(-1) ] ;
alpha=0;
beta=N;
ll = 15;
L = (N+1)*ll+1;

N_xmk = b_spline_evaluate(x_o, x_ext, alpha, beta, m, N, L);
dN_xmk_1 = db_spline_evaluate(1, x_o, x_ext, alpha, beta, m, N, L);
dN_xmk_2 = db_spline_evaluate(2, x_o, x_ext, alpha, beta, m, N, L);
dN_xmk_3 = db_spline_evaluate(3, x_o, x_ext, alpha, beta, m, N, L);

bd_spline_3 = [N_xmk(1,1), dN_xmk_1(1,1), dN_xmk_2(1,1), dN_xmk_3(1,1)];
bd_spline_2 = [N_xmk(1,2), dN_xmk_1(1,2), dN_xmk_2(1,2), dN_xmk_3(1,2)];
bd_spline_1 = [N_xmk(1,3), dN_xmk_1(1,3), dN_xmk_2(1,3), dN_xmk_3(1,3)];
int_spline_0 = [N_xmk(1,4), dN_xmk_1(1,4), dN_xmk_2(1,4), dN_xmk_3(1,4)];
A = [bd_spline_3; bd_spline_2; bd_spline_1; int_spline_0];

% create table of derivative values
fid = fopen('spline_derivatives_m4.txt', 'w');
% print a title, followed by a blank line
fprintf(fid, 'LHS boundary spline derivatives -- m=4 \n\n');
fprintf(fid, 'Rows: derivatives from order 0 -- 3 \n\n');
% print column headings
fprintf(fid, '%s %s %s %s\n', 'N_{-3}', 'N_{-2}', 'N_{-1}', 'N_{0}');
% print values in column order; two values appear on each row of the file
fprintf(fid, '%f  %f %f %f\n', A);
fclose(fid);