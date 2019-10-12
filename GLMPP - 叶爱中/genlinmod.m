function [y,yhat] = genlinmod(a,b,x,n1,n2,nmem);
% COMPUTE:
%    YHAT = AX
%    Y = YHAT + B*E, WHERE E IS A VECTOR OF SND'S
%       x is a column vector with n2 rows
%       a is an n1 x n2 array of coefficients
%
    [ma,na] = size(a);
    [mb,nb] = size(b);
    if (n1>0);
        if (n2>0);
           if (ma~=n1|na~=n2|mb~=n1|nb~=n1);
               yhat = 0;
               y = 0;
           else;
               yhat = a*x;
               for k=1:nmem;
                   e = normrnd(0,1,n1,1);
                   y(:,k) = yhat + b*e;
               end;
           end;
       else;
           yhat = zeros(n1:1);
           for k=1:nmem;
               e = normrnd(0,1,n1,1);
               y(:,k) = b*e;
           end;
       end;
    else;
        y = 0;
        yhat = 0;
    end;