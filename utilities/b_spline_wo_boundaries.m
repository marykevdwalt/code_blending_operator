function [Nxmkt] = b_spline_wo_boundaries(x, m, k, t)

%% joe's code
% 	oldies = zeros(1,m);
% 	for kk=1:m
% 		oldies(kk) = (t < x(kk-1+2)) & (t >= x(kk-1+1));
% 	end
% 		
% 	for mm=2:m
% 		newbies = zeros(1,m+1-mm);
% 		for kk=1:(m+1-mm)
% 			newbies(kk) = (t - x(kk))/(x(kk-1+mm) - x(kk))*oldies(kk) + (x(kk+mm) - t)/(x(kk+mm) - x(kk+1))*oldies(kk+1);
% 		end
% 		oldies = newbies;
% 	end
% 	Nxmkt = oldies(1);

%% adapted
if m==1
    if t < x(k+1)
        Nxmkt = 0;
    elseif t > x(k+2)
        Nxmkt = 0;
    else
        Nxmkt = 1;
    end
else
    TDD = divdiff(x(k+1:k+m+1), max(0, x(k+1:k+m+1) - t).^(m-1));
    Nxmkt = (x(k + m +1) - x(k+1)) * TDD(1, end) ;
end

end

	
