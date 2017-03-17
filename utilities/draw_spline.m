function [] = draw_spline(x_o,m,L,i)

x_ext=[ x_o(1) + [-m+1:-1]*1e-3 , x_o , x_o(end) + [1:m-1]*1e-3 ] ;
N=length(x_o)-2;
dx = linspace(x_o(i+1), x_o(i+m+1), L);
N_i = zeros(L,1);

for j = 1: L
	N_i(j) = b_spline(x_ext, m, N, i, dx(j)) ;
end

% figure
% plot(dx,N_i,'k', 'linewidth', 1.5);
% XTickLabel = {'$x_i$', '$x_{i+1}$', '$x_{i+2}$', '$x_{i+3}$', '$x_{i+4}$'};
% set(gca,'XTick',[x_o(i+1):1:x_o(i+1+m)]);
% set(gca,'YTick',[]);
% [hx,hy] = format_ticks(gca,XTickLabel, [],[],[],[],[],[],'FontSize', 40);

%fprintf('x(i-1) = %f \n', x_o(i-1));
%fprintf('x(i+3+m) = %f \n', x_o(i+3+m));

figure
plot(dx,N_i,'k', 'linewidth', 1.5); hold on;
plot([x_o(i : 2 : i+2+m)], zeros(size(1,2+m/2)), 'ko', 'markersize', 20);
xlim([x_o(i-1)-0.5 x_o(i+3+m)+0.5]);
set(gca,'XTick',[x_o(i-1 : 2 : i+3+m)]);
set(gca, 'TickLength', [0.025 0.025]);
set(gca,'YTick',[]);
XTickLabel = {'$\tilde{x}_{2i-2}$', '$\tilde{x}_{2i}$', '$\tilde{x}_{2i+2}=x_{i}$', '$\tilde{x}_{2i+4}$', '$\tilde{x}_{2i+6}$'};
[hx,hy] = format_ticks(gca,XTickLabel);

end