function [xmesh,cdf,z] = empcdf2(x);
% GENERATE EMPIRICAL CDF OF X
    n = length(x);
    [xsort,idx] = sort(x);
    for i=1:n-1;
        if xsort(i+1)<=xsort(i);
            xsort(i+1) = xsort(i) + .00001;
        end;
    end;
    range = xsort(n) - xsort(1);
    xmin = xsort(1) - range/10;
    xmax = xsort(n) + range/10;
    ncdf = n + 2;
    xmesh = zeros(ncdf,1);
    cdf = zeros(ncdf,1);
    xmesh(1) = xmin;
    xmesh(ncdf) = xmax;
    cdf(1) = 0.;
    cdf(ncdf) = 1.;
    for i=1:n;
      pr = (i-0.5)/n;
      cdf(i+1) = pr;
      xmesh(i+1) = xsort(i);
      z(idx(i),1) = norminv(pr);   %将经验分布转成标准正太分布
    end;
    n = n;
    