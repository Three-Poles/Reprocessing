function [qobs,qsim,nens,qsim2_ensm,a,b] = enspost_sim (qobs_calb,qsim_calb,nts,nts1,na,nf,ndays,jday,jmo,nmem,buffer);
%  ENSPOST_SIM = ROUTINE TO GENERATE ENSEMBLES OF QSIMX FOR SET OF QSIMX;
% Firstly use the first nts year to calibrate the GLMPP model parameters:
%1.  SET-UP ARRAYS QOBS AND QSIM;
%2.  MAP QOBS AND QSIM TO ZOBS, ZSIM AND ZERR by normal quantile transform 
%3. Fit GLMPP model
%Then use (nts+nts1) years to do post-processing and generate post-processed ensemble forecasts (after line 125)
%1.  SET-UP ARRAYS QOBS AND QSIM;
%2.  MAP QOBS AND QSIM TO ZOBS, ZSIM AND ZERR;
%3. Generate ensemble forecasts
% Author: Aizhong Ye, azye@bnu.edu.cn
 

nens = 0;
nbuff = buffer;

[ce,mobs,nobs] = size(qobs_calb);
[ce,msim,nsim] = size(qsim_calb);
if (mobs~=msim); return; end;
if (nobs~=nsim); return; end;
if (mobs<nts); return; end;
if (nobs<ndays); return; end;
%1.  SET-UP ARRAYS QOBS AND QSIM;
ndp = na + nf;
for i=1:na;  % 所给数据天数
    nens = 0;
    for k=1:nts1;  % 年数
        
        for j=1:nbuff;
            nens = nens + 1;  %年数×缓冲天数
            qobs(nens,i) = qobs_calb(1,k,i+j-1); %（每天采用的数据，天）
            qsim(nens,i) = qsim_calb(1,k,i+j-1); %585*24  39*39
        end;
    end;
end;
for i=na+1:ndp;  % 所给数据天数
    nens = 0; ce=i-na+1;
    for k=1:nts1;  % 年数
        njj=max(nbuff,1);
        for j=1:njj;
            nens = nens + 1;  %年数×缓冲天数
            qobs(nens,i) = qobs_calb(ce,k,i+j-1); %（每天采用的数据，天）
            qsim(nens,i) = qsim_calb(ce,k,i+j-1); %585*24  39*39
        end;
    end;
end;



%2.  MAP QOBS AND QSIM TO ZOBS, ZSIM AND ZERR; (normal quantile transform)
for i=1:ndp;
    [zsim(:,i),xmesh_qsimx(i,:),cdf_qsim(i,:),qsimx(:,i),bexp] = q2z(qsim(:,i));
    [zobs(:,i),xmesh_qobsx(i,:),cdf_qobs(i,:),qobsx(:,i),bexp] = q2z(qobs(:,i));
end;
zerr = zobs - zsim;
%  SET-UP DATA VECTORS
ndp = na + nf;
zsim1 = zsim(:,na+1:ndp);
zobs1 = zobs(:,na+1:ndp);
qsim1 = qsim(:,na+1:ndp);
qobs1 = qobs(:,na+1:ndp);
zerr1 = zobs1 - zsim1;
qerr1 = qsim1 - qobs1;
if (na>0);
    zsim2 = zsim(:,1:na); %（每天采用的数据，天）
    zobs2 = zobs(:,1:na);
    qsim2 = qsim(:,1:na);
    qobs2 = qobs(:,1:na);
else;
    zsim2 = 0.;
    zobs2 = 0.;
    qsim2 = 0.;
    qobs2 = 0.;
end;
corr_zobs = corr(zobs1); %观测自相关系数
corr_zsim = corr(zsim1);  % 模拟自相关系数
cross_corr_rawq = corr(qsim1,qobs1);
cross_corr_rawz = corr(zsim1,zobs1);
for j=1:nf;
    rmse_rawq(j) = rmserror(qsim1(:,j),qobs1(:,j));
end;
%
corr_zerr = corr(zerr1);

% GET STATISTICS OF QOBS, QSIM, ZOBS, ZSIM, ZERR, ZSIM_ENS AND QSIM_ENS
qsim_avg = mean(qsim1);
qsim_std = std(qsim1);
qobs_avg = mean(qobs1);
qobs_std = std(qobs1);
zsim_avg = mean(zsim1);
zsim_std = std(zsim1);
zobs_avg = mean(zobs1);
zobs_std = std(zobs1);
zerr_avg = mean(zerr1);
zerr_std = std(zerr1);
% GET RMSE OF QSIM
for j=1:nf;
    rmse_qsim(j) = rmserror(qsim1(:,j),qobs1(:,j));
    rmse_zsim(j) = rmserror(zsim1(:,j),zobs1(:,j));
    mae_qsim(j) = mae(qsim1(:,j),qobs1(:,j));
    mae_zsim(j) = mae(zsim1(:,j),zobs1(:,j));
end;
% SET UP ERROR MODEL 2
if (na>0);
    n = 2*(na + nf);
    z = zeros(nens,n);
    z = [zobs1 zsim1 zobs2 zsim2];
else
    n = 2*nf;
    z = zeros(nens,n);
    z = [zobs1 zsim1];
end;
clear a;
clear b;
clear bbt;
clear bbt_err;
n1 = nf;
n2 = nf + 2*na;
%3. fit GLMPP model
[a,b,bbt,bbt_err,d] = genlinmod_calb (n1,n2,z);%fit GLMPP model




%剔除缓冲区，重建Z进行后预报计算
%1.  SET-UP ARRAYS QOBS AND QSIM;
%clear qobs;
%clear qsim;
ndp = na + nf;
for i=1:na;  % 所给数据天数
    nens = 0;
    for k=1:nts;  % 年数
        for j=floor(nbuff/2)+1:floor(nbuff/2)+1;
            nens = nens + 1;  %年数×缓冲天数
            qobsf(nens,i) = qobs_calb(1,k,i+j-1); %（每天采用的数据，天）
            qsimf(nens,i) = qsim_calb(1,k,i+j-1); %585*24  39*39
        end;
    end;
end;

for i=na+1:ndp;  % 所给数据天数
    nens = 0;  kk=i-na+1;
    for k=1:nts;  % 年数
        for j=floor(nbuff/2)+1:floor(nbuff/2)+1;
            nens = nens + 1;  %年数×缓冲天数
            qobsf(nens,i) = qobs_calb(kk,k,i+j-1); %（每天采用的数据，天）
            qsimf(nens,i) = qsim_calb(kk,k,i+j-1); %585*24  39*39
        end;
    end;
end;


qobsf=[qobsf;qobs];%combine the two parts of data
qsimf=[qsimf;qsim];



%2.  MAP QOBS AND QSIM TO ZOBS, ZSIM AND ZERR;
clear zsim;
clear zobs;

for i=1:ndp;
    [zsim(:,i),xmesh_qsimx1(i,:),cdf_qsim1(i,:),qsimx1(:,i),bexp] = q2z(qsimf(:,i));
    [zobs(:,i),xmesh_qobsx1(i,:),cdf_qobs1(i,:),qobsx1(:,i),bexp] = q2z(qobsf(:,i));
end;
zerr = zobs - zsim;
%  SET-UP DATA VECTORS
ndp = na + nf;
zsim1 = zsim(:,na+1:ndp);
zobs1 = zobs(:,na+1:ndp);
qsim1 = qsimf(:,na+1:ndp);
qobs1 = qobsf(:,na+1:ndp);
zerr1 = zobs1 - zsim1;
qerr1 = qsim1 - qobs1;
if (na>0);
    zsim2 = zsim(:,1:na); %（每天采用的数据，天）
    zobs2 = zobs(:,1:na);
    qsim2 = qsimf(:,1:na);
    qobs2 = qobsf(:,1:na);
else;
    zsim2 = 0.;
    zobs2 = 0.;
    qsim2 = 0.;
    qobs2 = 0.;
end;
if (na>0);
    n = 2*(na + nf);
    z = zeros(nens,n);
    z = [zobs1 zsim1 zobs2 zsim2];
else
    n = 2*nf;
    z = zeros(nens,n);
    z = [zobs1 zsim1];
end;


%3. Generate zsim2_ens;
zsim2_ensm = zeros(nens,nf,nmem);
n = 2*(na + nf);
for i=1:nens;
    clear x;
    clear y;
    clear yhat;
    x(:,1) = z(i,nf+1:n);
    [y,yhat] = genlinmod(a,b,x,n1,n2,nmem);
    zsim2_ensm(i,:,:) = y;
end;
% Generate qsim2_ens
qsim2_ensm = zeros(nens,nf,nmem);

t = 1:nf;
nn = nens*nf;
cdfest = zeros(nn,1);

for i=1:nens;
    for k=1:nf;
        for j=1:nmem;
            %if (zsim2_ensm(i,k,j)<0) zsim2_ensm(i,k,j)=0; end;
            cdf2_ens(i,k,j) = normcdf(zsim2_ensm(i,k,j),0,1);
            qsimx2_ensm(i,k,j) = interp1q(cdf_qobs1(k+na,:)',xmesh_qobsx1(k+na,:)',cdf2_ens(i,k,j));
            qsim2_ensm(i,k,j) = qsimx2_ensm(i,k,j).^(1/bexp);
        end;
    end;
end;

n = 0;

