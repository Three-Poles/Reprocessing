function [rmse] = rmserror(x,y);
% COMPUTE RMSERROR BETWEEN X AND Y
    nx = length(x);
    ny = length(y);
    if (nx~=ny);
        rmse = -1;
    else;
        rmse = 0;
        for i=1:nx;
            rmse = rmse + (x(i)-y(i))^2;
        end;
        rmse = sqrt(rmse/nx);
    end
    
    