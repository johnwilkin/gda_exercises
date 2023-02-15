% For the CACSST example, examine the varimax rotated EOFs.
%
% The EOF analysis of CAC (Reynolds) SST data has to have been done by
% running script cacsst.m

colsst = cmap_anomghyde(128,0.02);
colsst = cmap_cmocean('curl');

if exist('D','var')~=1
  disp('Running script demo_1_cac_sst_eofs to do the EOF analysis')
  pauseit = 0;
  close all
  demo_1_cac_sst_eofs
  close all
end
close all
pauseit = 2;

% Retain only the signifcant EOFs because it makes no sense to try to apply
% physical interpretation to modes we have determined have no statistical
% significance.
K = 1:4;
Ut = U(:,K);
A = S*V';
At = A(K,:);

% Approximate data reconstruction for K modes is:
Dt = Ut*At;

% Notice that the Ut and At show that the EOF loadings are orthogonal
disp('EOF loading pattern correlations')
Ut'*Ut
disp(' ')

% ... and the amplitude time series are uncorrelated:
disp('EOF time series correlations')
At*At'

% compare to S^2
disp('S^2')
S(K,K).^2

%% Compute the varimax rotation of Ut
[Ur,R] = varimax_kaplan(Ut);

% The rotated matrix is Ur = Ut*R
% Check:
range(Ur-Ut*R)

% The new loadings (the columns of Un) are still orthogonal. Check:
disp('The new loadings (the columns of Un) are still orthogonal. Check:')
Ur'*Ur

% The (filtered) data matrix expansion can be written:
% Dt = Ut*R*inv(R)*At
% so we can write this as Dt = Ur*Ar with Ar = inv(R)*At;
Ar = inv(R)*At;

% Compare the loading patterns
for m=1:max(K)
  figure
  colormap(colsst)
  subplot(211)
  pcolorjw(lon,lat,reshape(U(:,m),[nlats nlons]))
  caxiscenter
  cax = caxis; colorbar
  title(['Original mode ' int2str(m)])
  subplot(212)
  pcolorjw(lon,lat,reshape(Ur(:,m),[nlats nlons]))
  caxis(cax); colorbar
  title(['Varimax rotated mode ' int2str(m)])
  pause_it
end

%% Compare %-variance explained for each mode
for m=1:3
  Dt = Ut(:,m)*At(m,:);
  Dr = Ur(:,m)*Ar(m,:);
  tmp = Dt-D;
  skt = 1-diag(tmp*tmp')./diag(D*D');
  tmp = Dr-D;
  skr = 1-diag(tmp*tmp')./diag(D*D');
  figure
  colormap(cmap_cmocean('tempo'))
  hax(1)=subplot(211);
  pcolorjw(lon,lat,reshape(100*skt,[nlats nlons]));
  cax = caxis; colorbar
  title(['% variance explained by mode ' int2str(m)])
  hax(2)=subplot(212);
  pcolorjw(lon,lat,reshape(100*skr,[nlats nlons]));
  caxis(cax); colorbar
  title(['% variance explained by rotated mode ' int2str(m)])
  linkaxes(hax);
  pause_it
end

%% Compare the amplitude time series

N = length(year);
for m=1:max(K)
  figure
  plot(year,At(m,:)/sqrt(N),year,Ar(m,:)/sqrt(N))
  legend(['Original mode ' int2str(m)],...
    ['Varimax rotated mode ' int2str(m)])
  grid on
  pause_it
end

% The total variance of the data is retained, but is distributed differently
% among the modes and has lost the property of maximizing the variance
% explained in the lowest mode, i.e. variance of rotated mode 1 is always
% less than unrotated.
tS = trace(S*S');
diag(At*At')/tS
diag(Ar*Ar')/tS
pause_it

%% And now the amplitude time series are correlated
% Compare:
At*At'
Ar*At'
pause_it

figure
subplot(211)
plot(year,At(1,:)/sqrt(N),year,At(3,:)/sqrt(N))
title('Modes 1 and 3 original')
grid on
subplot(212)
plot(year,Ar(1,:)/sqrt(N),year,Ar(3,:)/sqrt(N))
title('Modes 1 and 3 rotated')
grid on
pause_it

figure
[xct,~] = xcov(At(1,:),At(3,:));
[xcr,lags] = xcov(Ar(1,:),Ar(3,:));
plot(lags/12,xct,lags/12,xcr)
grid on
legend('x-cov modes 1 and 3 ','x-cov modes 1 and 3 rotated')
pause_it
disp('peak shifted toward left implies mode 1 leads mode 3')

%%
figure
colormap(colsst)
subplot(321)
m=1;
pcolorjw(lon,lat,reshape(Ut(:,m),[nlats nlons]))
caxiscenter
cax = caxis; colorbar
title(['Original mode ' int2str(m)])

subplot(323)
m=1;
pcolorjw(lon,lat,reshape(Ur(:,m),[nlats nlons]))
caxis(cax); colorbar
title(['Rotated mode ' int2str(m)])

subplot(322)
m=2;
pcolorjw(lon,lat,reshape(Ut(:,m),[nlats nlons]))
caxiscenter
cax = caxis; colorbar
title(['Original mode ' int2str(m)])

subplot(324)
m=2;
pcolorjw(lon,lat,reshape(Ur(:,m),[nlats nlons]))
caxis(cax); colorbar
title(['Rotated mode ' int2str(m)])

subplot(325)
[xcr,~] = xcov(Ar(1,:),Ar(2,:));
[xct,lags] = xcov(At(1,:),At(2,:));
plot(lags/12,xct,lags/12,xcr)
set(gca,'xlim',[-15 15])
grid on
legend('x-cov modes 1 and 2 ','x-cov modes 1 and 2 rotated')




