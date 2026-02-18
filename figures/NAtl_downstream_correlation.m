%% add all necessary paths

addpath /atalanta/home/swijffels/work/argo/
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/toolbox/cmocean
addpath /atalanta/home/chersh/toolbox
addpath /atalanta/home/swijffels/work/argo/gridNSF/strmat_sla_200rmax
addpath /atalanta/home/chersh/toolbox/cmocean
cd /atalanta/home/swijffels/work/argo/gridNSF

%{
addpath /home/swijffels/work/argo/
addpath /home/swijffels/toolbox/seawater
addpath /home/swijffels/toolbox/csirolib
addpath /home/swijffels/work/argo/matlab
addpath /home/swijffels/toolbox/subaxis
addpath /home/chersh/toolbox/cmocean
addpath /home/chersh/toolbox
addpath /home/swijffels/work/argo/gridNSF/strmat_sla_200rmax
addpath /home/chersh/toolbox/cmocean
cd /home/swijffels/work/argo/gridNSF
%}

%% load in argo climatology
if ~exist('clim')
    clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end
    %clim=load('gridSonSig1grid_climatology_fine_huber_sla_bar_xseas_1000.mat'); %load('gridonSig1grid_climatology_fine_huber_1000_sla.mat');
a_mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,1:85)) + squeeze(clim.pr.c(1:4:end,21:4:621,1:85,1));

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/atalanta/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end

%% get ecco data file
e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');   
%variables are: dnum, pri (pressure), pti (potential temperature), readme, si (salinity-181x141x241x312) , sig1grid, 
%ui, vi, xi, yi, yrgrid
%ep = matfile('/atalanta/home/swijffels/work/argo/gridNSF/ecco4/gridSonPgrid_ecco4r4_2020.mat');
e_sig1grid = e.sig1grid; 
e_yrgrid = e.yrgrid;
e_xi = e.xi;
e_yi = e.yi;

%% get argo data file
a = matfile('/atalanta/home/chersh/SpiceAnomalies/despike_argo_data.mat');
a_sig1grid = a.sig1grid; %sigma1 level grid
a_yrgrid = a.yrgrid; %months
a_xi = a.xi; %longitude grid
a_yi = a.yi; %latitude grid

%% pick a sigma1 level, load in argo and ECCO data slice
sigma1 =  '30.60'; %WHICH SIGMA LEVEL DO YOU WANT TO LOOK AT?
[~,e_isig]=min(abs(e_sig1grid - str2double(sigma1)));
[~,a_isig]=min(abs(a_sig1grid - str2double(sigma1))); %a_sig1grid and e_sig1grid aren't the same length

ii = 142; %WHICH STREAMLINE DO YOU WANT TO LOOK AT?
maxlag = 96; %months

a_sa = squeeze(a.sa_fwa(:,:,a_isig,:)); %argo salinity
a_dpr = squeeze(a.dpr_fwa(:,:,a_isig,:)); %argo layer thickness

e_si = squeeze(e.si(:,:,e_isig,:)); %pull out ecco salinity on this sigma1 surface
e_drhodri = squeeze(e.drhodri(:,:,e_isig,:)); %pull out ecco drho/dr on this sigma1 surface

load(['ECCO4r4_strm_ventilated_',num2str(sigma1),'.mat']); %load in mean ecco streamlines on this surface

%% find downstream correlations in Argo and ECCO

%SPICE

lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %longitude values of streamline
lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat)); %latitude values of streamline
cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %cumulative distance in km

%get Argo data along streamline
a_spice_mean = squeeze(clim.sa.m(:,:,a_isig)) + squeeze(clim.sa.c(:,:,a_isig,1)); %map of mean spice on this sigma surface
a_spice_mean = a_spice_mean(1:4:end,21:4:621); %climatology is on a finer grid than the anomalies

a_yrgrid_trunc = a_yrgrid(a_yrgrid >= 2004); %cut off before 2004
a_strmln_spice_an = NaN(length(lon_nums),length(a_yrgrid_trunc)); %spice/salinity anomalies along streamline
a_strmln_spice = NaN(length(lon_nums),length(a_yrgrid_trunc)); %absolute spice/salinity along streamline
a_strmln_spice_mean  = interp2(a_xi,a_yi,squeeze(a_spice_mean)',lon_nums,lat_nums,'linear'); %mean salinity along streamline

for it = 1:length(a_yrgrid_trunc)
   % fill out array of spice and spice anomaly from mean, along ecco streamline
   a_strmln_spice_an(:,it) = interp2(a_xi,a_yi,squeeze(a_sa(:,:,it+12))',lon_nums,lat_nums,'linear')-a_strmln_spice_mean;
   a_strmln_spice(:,it) = interp2(a_xi,a_yi,squeeze(a_sa(:,:,it+12))',lon_nums,lat_nums,'linear');
end

% ARGO: find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline
%griddata

%find the first point along the argo streamline that has no NaNs, then go a
%little past that
sp = 1;
while any(isnan(a_strmln_spice(sp,:)))
    sp = sp + 1;
end

startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
starting_spice = a_strmln_spice(startingpoint,:); %or strmln_spice_an
a_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
a_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
a_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location -- ensure these are negative so line doesn't jump around
a_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
for pp = 1:length(lon_nums) %for each point along the streamline
    downstream_spice = a_strmln_spice(pp,:); %or strmln_spice_an (anomaly from overall mean)
    if any(isnan(downstream_spice))
        continue
    end
    [r,lags] = xcorr(starting_spice,downstream_spice,maxlag,'coeff');
    a_R(:,pp) = r;
    a_Lags(:,pp) = lags;
    a_corrpeak(pp) = max(r(lags <= 0)); %% make sure lag is negative
    if pp > 100
        a_corrpeak(pp) = max(r(lags <= -20)); %% make sure lag is negative
    end
    ia = find(r == a_corrpeak(pp));
    a_advect(pp) = lags(ia);
    if sum(~isnan(r)) == 0
        a_advect(pp) = NaN;
    end
end

%calculate apparent speed of anomaly propagation, plus error bars (standard
%deviation of max correlation? )
%convert lag of max correlation to advection speed
%or use mean advective speed 
days_nums = ostrme(ii).days(~isnan(ostrme(ii).days)); %ecco streamline days since start, based on mean advective speed
%for some reason ecco days not always same length as other variables...
if length(days_nums) < length(cdxkm_nums)
    %lastdays = days_nums(end)-days_nums(end-1); %add one more day entry with same spacing as previous time step
    days_nums = [days_nums;zeros((length(cdxkm_nums)-length(days_nums)),1)];
end

months_nums = (days_nums-days_nums(startingpoint))/30.42;

a_meanadvectpeak = NaN(1,length(lon_nums));
a_midx = NaN(1,length(lon_nums));
for pp = 1:length(lon_nums)
    diff = NaN(length(lags),1);
    for ll = 1:length(lags)
        diff(ll) = abs(lags(ll) - months_nums(pp));
    end
    [~,a_midx(pp)] = min(diff);
    a_meanadvectpeak(pp) = a_R(length(lags)-a_midx(pp)+1,pp);
end

% ECCO: find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline
e_yrgrid_trunc = e_yrgrid(e_yrgrid >= 1993);
e_strmln_spice_an = NaN(length(lon_nums),length(e_yrgrid_trunc)); %spice/salinity anomalies along streamline
e_strmln_spice = NaN(length(lon_nums),length(e_yrgrid_trunc)); %absolute spice/salinity along streamline
e_spice_mean = nanmean(e_si,3); %map of mean salinity on this density surface
e_strmln_spice_mean  = interp2(e_xi,e_yi,squeeze(e_spice_mean)',lon_nums,lat_nums,'linear'); %mean salinity along streamline

for it = 1:length(e_yrgrid_trunc)
   % fill out array of spice and spice anomaly from mean, along ecco streamline
   e_strmln_spice_an(:,it) = interp2(e_xi,e_yi,squeeze(e_si(:,:,it))',lon_nums,lat_nums,'linear')-e_strmln_spice_mean;
   e_strmln_spice(:,it) = interp2(e_xi,e_yi,squeeze(e_si(:,:,it))',lon_nums,lat_nums,'linear');
end

%for each point along streamline, remove seasonal cycle
e_strmln_spice_rm_seasonal = NaN(length(lon_nums),length(e_yrgrid_trunc)); %spice/salinity anomalies along streamline
for is = 1:size(e_strmln_spice,1)
    if sum(~isnan(e_strmln_spice(is,:))) == 0 %seasonal.m doesn't like when it's all NaN!
        continue
    end
    [e_strmln_spice_rm_seasonal(is,:),zmean,coef,zm] = seasonal(e_yrgrid_trunc,e_strmln_spice(is,:)');
end

% find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline

e_starting_spice = e_strmln_spice_rm_seasonal(startingpoint,:); %use same startingpoint as for argo
e_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
e_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
e_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
e_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
for pp = 1:length(lon_nums) %for each point along the streamline
    e_downstream_spice = e_strmln_spice_rm_seasonal(pp,:); %or strmln_spice_an (anomaly from overall mean)
    [r,lags] = xcorr(e_starting_spice,e_downstream_spice,maxlag,'coeff');
    e_R(:,pp) = r;
    e_Lags(:,pp) = lags;
    e_corrpeak(pp) = max(r(lags <= 0)); %% make sure lag is negative
    if pp > 100
        e_corrpeak(pp) = max(r(lags <= -20)); %% make sure lag is negative
    end
    ia = find(r == e_corrpeak(pp));
    if size(ia) == [1,0]
        continue %weird, don't know why this is necessary
    end
    e_advect(pp) = lags(ia);
    
    if sum(~isnan(r)) == 0
        e_advect(pp) = NaN;
    end
end

%calculate apparent speed of anomaly propagation, plus error bars (standard
%deviation of max correlation? )
%convert lag of max correlation to advection speed
%or use mean advective speed from streamline files

e_meanadvectpeak = NaN(1,length(lon_nums));
e_midx = NaN(1,length(lon_nums));
for pp = 1:length(lon_nums)
    diff = NaN(length(lags),1);
    for ll = 1:length(lags)
        diff(ll) = abs(lags(ll) - months_nums(pp));
    end
    [~,e_midx(pp)] = min(diff);
    e_meanadvectpeak(pp) = e_R(length(lags)-e_midx(pp)+1,pp);
end

%plot mean salinity map with streamline overlayed
a_spice_mean(squeeze(a_mean_pr(:,:,a_isig)) <= clim_mldminmax.mldmax(1:4:end,21:4:621)) = NaN;
a_spice_mean(str2double(sigma1) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621)) = NaN;

figure('Position',[10 10 1000 700])
contourhighlights1 = [-0.5,0.5];
contourhighlights2 = [-0.7,0.7];
contourlevels = -1:0.1:1;
%figure(1)
%clf

subaxis(2,2,1,'SpacingVert',0.03,'SH',0.015) %make Argo correlation contour plot
%reverse the sign of lag on axis
contourf(-lags,cdxkm_nums,a_R',contourlevels,'edgecolor','none');
hold on
contour(-lags,cdxkm_nums,a_R',contourhighlights1,'k:','LineWidth',1);
contour(-lags,cdxkm_nums,a_R',contourhighlights2,'k','LineWidth',1);
plot([-maxlag maxlag],[cdxkm_nums(startingpoint) cdxkm_nums(startingpoint)],'r') %show which location is being correlated with
for ip = 1:floor(cdxkm_nums(end)/2000)
    [~,idx] = min(abs(cdxkm_nums - ip*2000));
    plot([-maxlag maxlag],[cdxkm_nums(idx) cdxkm_nums(idx)],'k:') %mark every 2000 km
end 
months_nums(283:end) = NaN;
a_advect(isnan(months_nums)) = NaN;
plot(-a_advect,cdxkm_nums,'k','LineWidth',2); %trace out lag of max correlation
plot(months_nums,cdxkm_nums,'b','LineWidth',2);
plot([0 0],[cdxkm_nums(1) cdxkm_nums(end)],'k');
ylabel('Downstream distance (km)')
set(gca,'XTick',[])
caxis([-1 1])
cmocean('balance',20)
vax = axis;
text(vax(1)+10,vax(4)-500,'Argo/spice')

subaxis(2,2,3) %make ECCO correlation contour plot
%contourf(lags,cdxkm_nums(1:(length(months_nums)/2)),e_R(:,1:(length(months_nums)/2))');
contourf(-lags,cdxkm_nums,e_R',contourlevels,'edgecolor','none');
hold on
contour(-lags,cdxkm_nums,e_R',contourhighlights1,'k:','LineWidth',1);
contour(-lags,cdxkm_nums,e_R',contourhighlights2,'k','LineWidth',1);
plot([-maxlag maxlag],[cdxkm_nums(startingpoint) cdxkm_nums(startingpoint)],'r') %show which location is being correlated with
for ip = 1:floor(cdxkm_nums(end)/2000)
    [~,idx] = min(abs(cdxkm_nums - ip*2000));
    plot([-maxlag maxlag],[cdxkm_nums(idx) cdxkm_nums(idx)],'k:') %mark every 2000 km
end
e_advect(isnan(months_nums)) = NaN;
plot(-e_advect,cdxkm_nums,'k','LineWidth',2); %trace out lag of max correlation
plot(months_nums,cdxkm_nums,'b','LineWidth',2);
plot([0 0],[cdxkm_nums(1) cdxkm_nums(end)],'k');
ylabel('Downstream distance (km)')
xlabel('Downstream lag (months)')
xticks([-80,-60,-40,-20,0,20,40,60,80])
caxis([-1 1])
cmocean('balance',20)
vax = axis;
text(vax(1)+10,vax(4)-500,'ECCO/spice');

%----------------------------------------------------------------------------------------------------------------
% PV: find downstream correlations in Argo and ECCO

rho_ref = 1000;

%get Argo data along this streamline
a_dpr_below = squeeze(a.dpr_fwa(:,:,a_isig+1,:));

%layer thickness climatology (pressure diff between sigma layer and the one above it)
a_dpclim2 = squeeze(clim.dpr.m(:,:,a_isig)) + squeeze(clim.dpr.c(:,:,a_isig,1));
%layer thickness climatology (pressure diff between sigma layer and the one below it)
a_dpclim1 = squeeze(clim.dpr.m(:,:,a_isig+1)) + squeeze(clim.dpr.c(:,:,a_isig+1,1));

%calculate fcor with clim dimensions
a_fcor = calc_fcor(clim.xi,clim.yi);

a_drho2 = a_sig1grid(a_isig) - a_sig1grid(a_isig-1); %change in density between this layer and the one above it
a_drho1 = a_sig1grid(a_isig+1) - a_sig1grid(a_isig); %change in density between this layer and the one below it

a_drhodpclim = (a_drho2./a_dpclim2 + a_drho1./a_dpclim1)/2;
a_pvclim = (1/rho_ref)*a_fcor.*a_drhodpclim; %pv climatological mean overall
a_pv =  NaN(size(a_dpr)); %calculate pv anomaly for argo

for it=1:length(a_yrgrid)
    a_dpr2 = a_dpr(:,:,it); %anomaly of layer thickness between this layer and one above it
    a_dpr1 = a_dpr_below(:,:,it); %anomaly of layer thickness between this layer and one below it
    a_drhodp_tot = (a_drho2./(a_dpclim2(1:4:end,21:4:621) + a_dpr2) + a_drho1./(a_dpclim1(1:4:end,21:4:621) + a_dpr1))/2;
    a_pv(:,:,it) = (1/rho_ref)*a_fcor(1:4:end,21:4:621).*a_drhodp_tot - a_pvclim(1:4:end,21:4:621);
end  

a_yrgrid_trunc = a_yrgrid(a_yrgrid >= 2004); %truncate before 2004
a_strmln_pv = NaN(length(lon_nums),length(a_yrgrid_trunc)); %pv anomaly along streamline

for it = 1:length(a_yrgrid_trunc)
   % fill out array of pv anomaly from mean, along ecco streamline
   a_strmln_pv(:,it) = interp2(a_xi,a_yi,squeeze(a_pv(:,:,it+12))',lon_nums,lat_nums,'linear'); %look into to see if somehow remaining nans where there shouldn't be
end

a_strmln_pv = a_strmln_pv(:,2:end); %for some reason the first month is NaNs

% ARGO: find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline
%find the first point along the streamline that has no NaNs, then go a
%little past that
%or, find the first point along the streamline that has 
sp = 1;
while any(isnan(a_strmln_pv(sp,:)))
    sp = sp + 1;
end
outcrop = sp;

startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
a_starting_pv = a_strmln_pv(startingpoint,:); %or strmln_spice_an
a_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
a_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
a_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
a_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
for pp = 1:length(lon_nums) %for each point along the streamline
    a_downstream_pv = a_strmln_pv(pp,:); %or strmln_spice_an (anomaly from overall mean)
    [r,lags] = xcorr(a_starting_pv,a_downstream_pv,maxlag,'coeff');
    a_R(:,pp) = r;
    a_Lags(:,pp) = lags;
    a_corrpeak(pp) = max(r(lags <= 0 & lags >= -60)); %% make sure lag is negative
    ia = find(r == a_corrpeak(pp));
    if size(ia) == [1,0]
        continue %weird, don't know why this is necessary
    end
    a_advect(pp) = lags(ia);
    if sum(~isnan(r)) == 0
        a_advect(pp) = NaN;
    end
end

%calculate apparent speed of anomaly propagation
%convert lag of max correlation to advection speed
%or use mean advective speed 
days_nums = ostrme(ii).days(~isnan(ostrme(ii).days)); %ecco streamline days since start, based on mean advective speed
%for some reason ecco days not always same length as other variables...
if length(days_nums) < length(cdxkm_nums)
    %lastdays = days_nums(end)-days_nums(end-1); %add one more day entry with same spacing as previous time step
    days_nums = [days_nums;zeros((length(cdxkm_nums)-length(days_nums)),1)];
end

months_nums = (days_nums-days_nums(startingpoint))/30.42;

a_meanadvectpeak = NaN(1,length(lon_nums));
a_midx = NaN(1,length(lon_nums));
for pp = 1:length(lon_nums)
    diff = NaN(length(lags),1);
    for ll = 1:length(lags)
        diff(ll) = abs(lags(ll) - months_nums(pp));
    end
    [~,a_midx(pp)] = min(diff);
    a_meanadvectpeak(pp) = a_R(length(lags)-a_midx(pp)+1,pp);
end

% ECCO: find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline
e_yrgrid_trunc = e_yrgrid(e_yrgrid >= 1993);
e_fcor = calc_fcor(e_xi,e_yi);
e_drhodrclim = nanmean(e_drhodri,3);
e_pvclim = (1/rho_ref)*e_fcor.*e_drhodrclim;
e_pv = -1*((1/rho_ref)*e_fcor.*e_drhodri - e_pvclim); %times -1 b/c drhodri is negative

e_strmln_pv = NaN(length(lon_nums),length(e_yrgrid_trunc)); %pv anomalies along streamline
e_strmln_pv_mean  = interp2(e_xi,e_yi,squeeze(e_pvclim)',lon_nums,lat_nums,'linear'); %mean pv along streamline

for it = 1:length(e_yrgrid_trunc)
   % fill out array of spice and spice anomaly from mean, along ecco streamline
   e_strmln_pv(:,it) = interp2(e_xi,e_yi,squeeze(e_pv(:,:,it))',lon_nums,lat_nums,'linear')-e_strmln_pv_mean;
end

%for each point along streamline, remove seasonal cycle
e_strmln_pv_rm_seasonal = NaN(length(lon_nums),length(e_yrgrid_trunc)); %spice/salinity anomalies along streamline
for is = 1:size(e_strmln_pv,1)
    if sum(~isnan(e_strmln_pv(is,:))) == 0 %seasonal.m doesn't like when it's all NaN!
        continue
    end
    [e_strmln_pv_rm_seasonal(is,:),zmean,coef,zm] = seasonal(e_yrgrid_trunc,e_strmln_pv(is,:)');
end

% find lagged correlation coefficients from point near outcrop (past
% seasonal cycle) and points farther down streamline

e_starting_pv = e_strmln_pv_rm_seasonal(startingpoint,:); %use same startingpoint as argo
e_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
e_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
e_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
e_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
for pp = 1:length(lon_nums) %for each point along the streamline
    e_downstream_pv = e_strmln_pv_rm_seasonal(pp,:); %or strmln_spice_an (anomaly from overall mean)
    [r,lags] = xcorr(e_starting_pv,e_downstream_pv,maxlag,'coeff');
    e_R(:,pp) = r;
    e_Lags(:,pp) = lags;
    e_corrpeak(pp) = max(r(lags <= 0 & lags >= -50)); %% make sure lag is negative
    ia = find(r == e_corrpeak(pp));
    if size(ia) == [1,0]
        continue %weird, don't know why this is necessary
    end
    e_advect(pp) = lags(ia);
end

%calculate apparent speed of anomaly propagation, plus error bars (standard
%deviation of max correlation? )
%convert lag of max correlation to advection speed
%or use mean advective speed from streamline files
e_meanadvectpeak = NaN(1,length(lon_nums));
e_midx = NaN(1,length(lon_nums));
for pp = 1:length(lon_nums)
    diff = NaN(length(lags),1);
    for ll = 1:length(lags)
        diff(ll) = abs(lags(ll) - months_nums(pp));
    end
    [~,e_midx(pp)] = min(diff);
    e_meanadvectpeak(pp) = e_R(length(lags)-e_midx(pp)+1,pp);
end

%plot mean pv map with streamline overlayed
%first, crop pvclim to winter outcrop
a_pvclim = a_pvclim(1:4:end,21:4:621);
a_pvclim(squeeze(a_mean_pr(:,:,a_isig)) <= clim_mldminmax.mldmax(1:4:end,21:4:621)) = NaN;
a_pvclim(str2double(sigma1) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621)) = NaN;

subaxis(2,2,2) %make argo correlation contour plot
contourf(-lags,cdxkm_nums,a_R',contourlevels,'edgecolor','none');
hold on
contour(-lags,cdxkm_nums,a_R',contourhighlights1,'k:','LineWidth',1);
contour(-lags,cdxkm_nums,a_R',contourhighlights2,'k','LineWidth',1);
plot([-maxlag maxlag],[cdxkm_nums(startingpoint) cdxkm_nums(startingpoint)],'r') %show which location is being correlated with
for ip = 1:floor(cdxkm_nums(end)/2000)
    [~,idx] = min(abs(cdxkm_nums - ip*2000));
    plot([-maxlag maxlag],[cdxkm_nums(idx) cdxkm_nums(idx)],'k:') %mark every 2000 km
end
months_nums(283:end) = NaN;
a_advect(isnan(months_nums)) = NaN;
plot(-a_advect,cdxkm_nums,'k','LineWidth',2); %trace out lag of max correlation
plot(months_nums,cdxkm_nums,'b','LineWidth',2);
plot([0 0],[cdxkm_nums(1) cdxkm_nums(end)],'k');
set(gca,'YTick',[])
set(gca,'XTick',[])
caxis([-1 1])
cmocean('balance',20)
vax = axis;
text(vax(1)+10,vax(4)-500,'Argo/PV')

subaxis(2,2,4) %make ecco correlation contour plot
contourf(-lags,cdxkm_nums,e_R',contourlevels,'edgecolor','none');
hold on
contour(-lags,cdxkm_nums,e_R',contourhighlights1,'k:','LineWidth',1);
contour(-lags,cdxkm_nums,e_R',contourhighlights2,'k','LineWidth',1);
plot([-maxlag maxlag],[cdxkm_nums(startingpoint) cdxkm_nums(startingpoint)],'r') %show which location is being correlated with
for ip = 1:floor(cdxkm_nums(end)/2000)
    [~,idx] = min(abs(cdxkm_nums - ip*2000));
    plot([-maxlag maxlag],[cdxkm_nums(idx) cdxkm_nums(idx)],'k:') %mark every 2000 km
end
e_advect(isnan(months_nums)) = NaN;
plot(-e_advect,cdxkm_nums,'k','LineWidth',2); %trace out lag of max correlation
plot(months_nums,cdxkm_nums,'b','LineWidth',2);
plot([0 0],[cdxkm_nums(1) cdxkm_nums(end)],'k');
xlabel('Downstream lag (months)')
xticks([-80,-60,-40,-20,0,20,40,60,80])
set(gca,'YTick',[])
c1 = colorbar;
c1.Label.String = 'correlation coefficient';
c1.Orientation = 'Horizontal';
c1.Position = [.3 .91 .4 .01];
caxis([-1 1])
cmocean('balance',20)
vax = axis;
text(vax(1)+10,vax(4)-500,'ECCO/PV')

saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_downstream_correlation.png');
saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_downstream_correlation.fig');

%% useful functions
function fcor = calc_fcor(xgrid,ygrid)
    fcor = NaN(length(xgrid),length(ygrid)); %coriolis parameter
    for aa = 1:length(xgrid)
        for bb = 1:length(ygrid)
            lati = ygrid(bb);
            fcor(aa,bb) = 2*7.2921*10^(-5)*sind(lati); %fill out coriolis parameter array
        end
    end
end