function [zsim,xmesh_qsimx,cdf_qsimx,qsimx,b] = q2z(qsim);
% MAP Q TO Z;
 
    b = sqrt(1/1);
    [m,n] = size(qsim);
    qsimx = qsim.^b;
    for j=1,n;
        [xmesh_qsimx(:,j),cdf_qsimx(:,j),zsim] = empcdf2(qsimx(:,j));
    end;
    n = n;

