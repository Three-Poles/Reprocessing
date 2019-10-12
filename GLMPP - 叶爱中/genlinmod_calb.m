function [a,b,bbt,bbt_err,d] = genlinmod_calb(n1,n2,z);
% CALIBRATE PARAMETERS OF A GENERAL LINEAR MODEL
%   Z is an input data array with standardized (0,1) observations 
%   Z has nobs rows and n1+n2 columns
%   The first n1 columns of z contain dependent variables
%   The last n2 columns of z contain independent variables
    [m,n] = size(z);
    if (n~=n1+n2);
        a = zeros(n1:1);
        bbt = zeros(n1:n1);
        bbt = zeros(n1:n1);
        bbt_err = zeros(n1:n1);
    elseif (n2==0);
        a = zeros(n1:1);
        bbt = cov(z,1);
    else;
        s = cov(z,1);
        n = n1 + n2;
        s11 = s(1:n1,1:n1);
        s12 = s(1:n1,n1+1:n);
        s21 = s(n1+1:n,1:n1);
        s22 = s(n1+1:n,n1+1:n);
        ta=cond(s22);
        if (ta< 1e+15); 
          s22inv = inv(s22);
          a = s12*s22inv;   
        else
          a = eye(n1,n-n1);  
          %a = zeros(n1,n-n1);  
        end
        
        bbt = s11 - a*s21;
    end;
    [v,d] = eig(bbt);
    dsqrt = sqrt(d);
    b = v*dsqrt;
    bbt_chk = b*b';
    bbt_err = bbt - bbt_chk;
    if (ta>1e+15); b = zeros(n1,n1);  end;
    %xxx=sum(sum(abs(bbt_err)))
    %if (xxx>1)
    %nnn=1;
    %end
