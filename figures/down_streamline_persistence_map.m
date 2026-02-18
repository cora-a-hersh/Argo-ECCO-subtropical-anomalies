%% make a map of persistence of spice anomalies down streamlines on an isopycnal

% add paths
addpath /atalanta/home/swijffels/work/argo/
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/toolbox/cmocean
addpath /atalanta/home/chersh/toolbox
addpath /atalanta/home/swijffels/work/argo/gridNSF/strmat_sla_200rmax
addpath /atalanta/home/chersh/toolbox/cmocean

%%notes
%also save cumulative distance for each point
%weed streamlines to look at just tunnel streamlines
%also save pressure of sigma surface at each point

% load in Argo climatology
if ~exist('clim')
    clim = load('/atalanta/home/swijffels/work/argo/gridNSF/sig1grid/gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end
    %clim=load('gridSonSig1grid_climatology_fine_huber_sla_bar_xseas_1000.mat'); %load('gridonSig1grid_climatology_fine_huber_1000_sla.mat');
a_mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,1:85)) + squeeze(clim.pr.c(1:4:end,21:4:621,1:85,1));

% get ecco & argo data files
e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');   
%variables are: dnum, pri (pressure), pti (potential temperature), readme, si (salinity-181x141x241x312) , sig1grid, 
%ui, vi, xi, yi, yrgrid
%ep = matfile('/atalanta/home/swijffels/work/argo/gridNSF/ecco4/gridSonPgrid_ecco4r4_2020.mat');
e_sig1grid = e.sig1grid; 
e_yrgrid = e.yrgrid;
e_xi = e.xi;
e_yi = e.yi;

a = matfile('/atalanta/home/chersh/SpiceAnomalies/despike_argo_data.mat');
a_sig1grid = a.sig1grid; %sigma1 level grid
a_yrgrid = a.yrgrid; %months
a_xi = a.xi; %longitude grid
a_yi = a.yi; %latitude grid

%% pick a sigma1 level, load in argo and ECCO data slice
sigma1 =  '29.55'; %WHICH SIGMA LEVEL DO YOU WANT TO LOOK AT?
[~,e_isig]=min(abs(e_sig1grid - str2double(sigma1)));
[~,a_isig]=min(abs(a_sig1grid - str2double(sigma1))); %a_sig1grid and e_sig1grid aren't the same length

a_sa = squeeze(a.sa_fwa(:,:,a_isig,:)); %argo salinity
a_dpr = squeeze(a.dpr_fwa(:,:,a_isig,:)); %argo layer thickness

e_si = squeeze(e.si(:,:,e_isig,:)); %pull out ecco salinity on this sigma1 surface
e_drhodri = squeeze(e.drhodri(:,:,e_isig,:)); %pull out ecco drho/dr on this sigma1 surface

load(['/atalanta/home/swijffels/work/argo/gridNSF/strmat_sla_200rmax/ECCO4r4_strm_ventilated_',num2str(sigma1),'.mat']); %load in mean ecco streamlines on this surface
sig_a_mean_pr = squeeze(a_mean_pr(:,:,a_isig));

%% ECCO and Argo: for each streamline, calculate persistence of spice and PV anomalies at each point
tunnels = readtable('/home/chersh/Downloads/tunnels.csv','NumHeaderlines',1);
%ECCO spice 
latpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of latitude coordinates
lonpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of longitude coordinates
distkm = NaN(length(ostrme),length(ostrme(1).lon)); %list of distances in km from starting point at surface outcrop
mean_pr = NaN(length(ostrme),length(ostrme(1).lon)); %list of mean pressure values at each point
max_corr = NaN(length(ostrme),length(ostrme(1).lon)); %max correlation at each location along streamline
advect_corr = NaN(length(ostrme),length(ostrme(1).lon)); %correlation following mean advective speed
max_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to max correlation at each location
advect_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to mean advective speed
starting_variance = NaN(length(ostrme),length(ostrme(1).lon)); %variance of salinity at beginning of streamline

maxlag = 96; %months
%for ii = 1:length(ostrme) %run through all streamlines on this surface
streamlines = tunnels.Var12; %Var# is index of the current sig1 surface in list of densities for which there are down streamline persistence files
streamlines = streamlines(~isnan(streamlines));
for ii = streamlines'
    disp(ii)
    if length(ostrme(ii).lat) ~= length(ostrme(1).lat)
        continue
    end
    latpoint(ii,:) = ostrme(ii).lat;
    lonpoint(ii,:) = ostrme(ii).lon;
    distkm(ii,:) = ostrme(ii).cdxkm;
    
    lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %longitude values of streamline
    lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat)); %latitude values of streamline
    if length(lat_nums) - length(lon_nums) == 1 %sometimes length of these are 1 off?
        lat_nums = lat_nums(1:(end-1));
    elseif length(lon_nums) - length(lat_nums) == 1
        lon_nums = lon_nums(1:(end-1));
    end
    cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %cumulative distance in km
    
    mean_pr(ii,1:length(lon_nums)) = interp2(a_xi,a_yi,sig_a_mean_pr',lon_nums,lat_nums,'linear');

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

    % ECCO: find lagged correlation coefficients from point near outcrop (past
    % seasonal cycle) and points farther down streamline
    e_strmln_spice_an = NaN(length(lon_nums),length(e_yrgrid)); %spice/salinity anomalies along streamline
    e_strmln_spice = NaN(length(lon_nums),length(e_yrgrid)); %absolute spice/salinity along streamline
    e_spice_mean = nanmean(e_si,3); %map of mean salinity on this density surface
    e_strmln_spice_mean  = interp2(e_xi,e_yi,squeeze(e_spice_mean)',lon_nums,lat_nums,'linear'); %mean salinity along streamline

    for it = 1:length(e_yrgrid)
    % fill out array of spice and spice anomaly from mean, along ecco streamline
        e_strmln_spice_an(:,it) = interp2(e_xi,e_yi,squeeze(e_si(:,:,it))',lon_nums,lat_nums,'linear')-e_strmln_spice_mean;
        e_strmln_spice(:,it) = interp2(e_xi,e_yi,squeeze(e_si(:,:,it))',lon_nums,lat_nums,'linear');
    end

    %for each point along streamline, remove seasonal cycle
    e_strmln_spice_rm_seasonal = NaN(length(lon_nums),length(e_yrgrid)); %spice/salinity anomalies along streamline
    for is = 1:size(e_strmln_spice,1)
        if sum(~isnan(e_strmln_spice(is,:))) < 5 %seasonal.m doesn't like when it's all NaN!
            continue
        end
        [e_strmln_spice_rm_seasonal(is,:),zmean,coef,zm] = seasonal(e_yrgrid,e_strmln_spice(is,:)');
    end

    %find the first point along the ecco streamline that has no NaNs, then go a
    %little past that
    sp = 1;
    %get an error when sp longer than 1st dim of e_strmln_spice_rm_seasonal
    stop = size(e_strmln_spice_rm_seasonal,1);
    while any(isnan(e_strmln_spice_rm_seasonal(sp,:))) && sp < stop
        sp = sp + 1;
    end
    if sp == stop
        continue
    end

    startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
   
    if startingpoint >= length(days_nums)
        continue
    else
        months_nums = (days_nums-days_nums(startingpoint))/30.42;
    end
    
    % find lagged correlation coefficients from point near outcrop (past
    % seasonal cycle) and points farther down streamline

    e_starting_spice = e_strmln_spice_rm_seasonal(startingpoint,:); %use same startingpoint as for argo
    starting_variance(ii,:) = var(e_starting_spice,'omitnan'); %calculate variance at beginning of streamline 
    e_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    e_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
    e_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
    for pp = 1:length(lon_nums) %for each point along the streamline
        e_downstream_spice = e_strmln_spice_rm_seasonal(pp,:); %or strmln_spice_an (anomaly from overall mean)
        [r,lags] = xcorr(e_starting_spice,e_downstream_spice,maxlag,'coeff');
        e_R(:,pp) = r;
        [e_corrpeak(pp),ia] = max(r(lags <= 0)); %%make sure lag is negative
        max_corr(ii,pp) = e_corrpeak(pp); %save correlations along max correlation line
        max_lag(ii,pp) = lags(ia);
        if sum(~isnan(r)) == 0
            max_lag(ii,pp) = NaN;
        end
    end

    %calculate apparent speed of anomaly propagation, plus error bars (standard
    %deviation of max correlation? )
    %convert lag of max correlation to advection speed
    %or use mean advective speed from streamline files
    
    e_midx = NaN(1,length(lon_nums));
    for pp = 1:length(lon_nums)
        diff = NaN(length(lags),1);
        for ll = 1:length(lags)
            diff(ll) = abs(lags(ll) - months_nums(pp));
        end
        [~,e_midx(pp)] = min(diff);
        advect_corr(ii,pp) = e_R(length(lags)-e_midx(pp)+1,pp);
        advect_lag(ii,pp) = lags(e_midx(pp)); 
    end
end

%turn (streamline,point) arrays into 1D vectors
latpoint = reshape(latpoint,[length(ostrme)*3000,1]);
lonpoint = reshape(lonpoint,[length(ostrme)*3000,1]);
distkm = reshape(distkm,[length(ostrme)*3000,1]);
mean_pr = reshape(mean_pr,[length(ostrme)*3000,1]);
max_lag =  reshape(max_lag,[length(ostrme)*3000,1]);
advect_lag = reshape(advect_lag,[length(ostrme)*3000,1]);
max_corr = reshape(max_corr,[length(ostrme)*3000,1]);
advect_corr = reshape(advect_corr,[length(ostrme)*3000,1]);
starting_variance = reshape(starting_variance,[length(ostrme)*3000,1]);

%remove points with no data!
keep = ~isnan(latpoint) & ~isnan(lonpoint);
latpoint = latpoint(keep);
lonpoint = lonpoint(keep);
distkm = distkm(keep);
mean_pr = mean_pr(keep);
max_lag =  max_lag(keep);
advect_lag = advect_lag(keep);
max_corr = max_corr(keep);
advect_corr = advect_corr(keep);
starting_variance = starting_variance(keep);

%save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/ecco_down_streamline_persistence_sig',sigma1,'.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');
save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/ecco_down_streamline_persistence_sig',sigma1,'_tunnels.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');

%-----------------------------------------------------------------------------------------------------
% ECCO: for each streamline, calculate persistence of PV anomalies at each point

latpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of latitude coordinates
lonpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of longitude coordinates
distkm = NaN(length(ostrme),length(ostrme(1).lon)); %list of distances in km from starting point at surface outcrop
mean_pr = NaN(length(ostrme),length(ostrme(1).lon)); %list of mean pressure values at each point
max_corr = NaN(length(ostrme),length(ostrme(1).lon)); %max correlation at each location along streamline
advect_corr = NaN(length(ostrme),length(ostrme(1).lon)); %correlation following mean advective speed
max_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to max correlation at each location
advect_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to mean advective speed
starting_variance = NaN(length(ostrme),length(ostrme(1).lon)); %variance of salinity at beginning of streamline

rho_ref = 1000;

maxlag = 96; %months
for ii = streamlines' %1:length(ostrme) %run through all streamlines on this surface
    disp(ii)
    if length(ostrme(ii).lat) ~= length(ostrme(1).lat)
        continue
    end
    latpoint(ii,:) = ostrme(ii).lat;
    lonpoint(ii,:) = ostrme(ii).lon;
    distkm(ii,:) = ostrme(ii).cdxkm;
    
    lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %longitude values of streamline
    lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat)); %latitude values of streamline
    if length(lat_nums) - length(lon_nums) == 1 %sometimes length of these are 1 off?
        lat_nums = lat_nums(1:(end-1));
    elseif length(lon_nums) - length(lat_nums) == 1
        lon_nums = lon_nums(1:(end-1));
    end
    cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %cumulative distance in km

    mean_pr(ii,1:length(lon_nums)) = interp2(a_xi,a_yi,sig_a_mean_pr',lon_nums,lat_nums,'linear');
    
    if length(lon_nums) < 20
       continue
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
    
    % ECCO: find lagged correlation coefficients from point near outcrop (past
    % seasonal cycle) and points farther down streamline
    %e_yrgrid_trunc = e_yrgrid(e_yrgrid >= 2004);
    e_fcor = calc_fcor(e_xi,e_yi);
    e_drhodrclim = nanmean(e_drhodri,3);
    e_pvclim = (1/rho_ref)*e_fcor.*e_drhodrclim;
    e_pv = -1*((1/rho_ref)*e_fcor.*e_drhodri - e_pvclim); %times -1 b/c drhodri is negative

    e_strmln_pv = NaN(length(lon_nums),length(e_yrgrid)); %pv anomalies along streamline
    e_strmln_pv_mean  = interp2(e_xi,e_yi,squeeze(e_pvclim)',lon_nums,lat_nums,'linear'); %mean pv along streamline

    for it = 1:length(e_yrgrid)
    % fill out array of spice and spice anomaly from mean, along ecco streamline
        e_strmln_pv(:,it) = interp2(e_xi,e_yi,squeeze(e_pv(:,:,it))',lon_nums,lat_nums,'linear')-e_strmln_pv_mean;
    end

    %for each point along streamline, remove seasonal cycle
    e_strmln_pv_rm_seasonal = NaN(length(lon_nums),length(e_yrgrid)); %spice/salinity anomalies along streamline
    for is = 1:size(e_strmln_pv,1)
        if sum(~isnan(e_strmln_pv(is,:))) < 6 %seasonal.m doesn't like when it's all NaN!
            continue
        end
        [e_strmln_pv_rm_seasonal(is,:),zmean,coef,zm] = seasonal(e_yrgrid,e_strmln_pv(is,:)');
    end

     %find the first point along the ecco streamline that has no NaNs, then go a
    %little past that
    sp = 1;
    %get an error when sp longer than 1st dim of e_strmln_spice_rm_seasonal
    stop = size(e_strmln_pv_rm_seasonal,1);
    while any(isnan(e_strmln_pv_rm_seasonal(sp,:))) && sp < stop
        sp = sp + 1;
    end
    if sp == stop
        continue
    end

    startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
   
    if startingpoint >= length(days_nums)
        continue
    else
        months_nums = (days_nums-days_nums(startingpoint))/30.42;
    end
    
    % find lagged correlation coefficients from point near outcrop (past
    % seasonal cycle) and points farther down streamline

    e_starting_pv = e_strmln_pv_rm_seasonal(startingpoint,:); %use same startingpoint as argo
    starting_variance(ii,:) = var(e_starting_pv,'omitnan'); %calculate variance at beginning of streamline 
    e_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    e_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    e_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
    e_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
    for pp = 1:length(lon_nums) %for each point along the streamline
        e_downstream_pv = e_strmln_pv_rm_seasonal(pp,:); %or strmln_spice_an (anomaly from overall mean)
        [r,lags] = xcorr(e_starting_pv,e_downstream_pv,maxlag,'coeff');
        e_R(:,pp) = r;
        e_Lags(:,pp) = lags;
        [max_corr(ii,pp),ia] =  max(r(lags <= 0)); %% make sure lag is negative
        max_lag(ii,pp) = lags(ia);
        if sum(~isnan(r)) == 0
            max_lag(ii,pp) = NaN;
        end
    end

    %use mean advective speed from streamline files
    e_midx = NaN(1,length(lon_nums));
    for pp = 1:length(lon_nums)
        diff = NaN(length(lags),1);
        for ll = 1:length(lags)
            diff(ll) = abs(lags(ll) - months_nums(pp));
        end
        [~,e_midx(pp)] = min(diff);
        advect_corr(ii,pp) = e_R(length(lags)-e_midx(pp)+1,pp);
        advect_lag(ii,pp) = lags(e_midx(pp)); 
    end
    
end

%turn (streamline,point) arrays into 1D vectors
latpoint = reshape(latpoint,[length(ostrme)*3000,1]);
lonpoint = reshape(lonpoint,[length(ostrme)*3000,1]);
distkm = reshape(distkm,[length(ostrme)*3000,1]);
mean_pr = reshape(mean_pr,[length(ostrme)*3000,1]);
max_lag =  reshape(max_lag,[length(ostrme)*3000,1]);
advect_lag = reshape(advect_lag,[length(ostrme)*3000,1]);
max_corr = reshape(max_corr,[length(ostrme)*3000,1]);
advect_corr = reshape(advect_corr,[length(ostrme)*3000,1]);
starting_variance = reshape(starting_variance,[length(ostrme)*3000,1]);

%remove points with no data!
keep = ~isnan(latpoint) & ~isnan(lonpoint);
latpoint = latpoint(keep);
lonpoint = lonpoint(keep);
distkm = distkm(keep);
mean_pr = mean_pr(keep);
max_lag =  max_lag(keep);
advect_lag = advect_lag(keep);
max_corr = max_corr(keep);
advect_corr = advect_corr(keep);
starting_variance = starting_variance(keep);

%save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/PV_ecco_down_streamline_persistence_sig',sigma1,'.mat'],'sigma1','latpoint','lonpoint','mean_pr','distkm','max_corr','advect_corr','max_lag','advect_lag','starting_variance');
save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/PV_ecco_down_streamline_persistence_sig',sigma1,'_tunnels.mat'],'sigma1','latpoint','lonpoint','mean_pr','distkm','max_corr','advect_corr','max_lag','advect_lag','starting_variance');

%--------------------------------------------------------------------------------------------------------
% Argo: for each streamline, calculate persistence of spice anomalies at each point

latpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of latitude coordinates
lonpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of longitude coordinates
distkm = NaN(length(ostrme),length(ostrme(1).lon)); %list of distances in km from starting point at surface outcrop
mean_pr = NaN(length(ostrme),length(ostrme(1).lon)); %list of mean pressure values at each point
max_corr = NaN(length(ostrme),length(ostrme(1).lon)); %max correlation at each location along streamline
advect_corr = NaN(length(ostrme),length(ostrme(1).lon)); %correlation following mean advective speed
max_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to max correlation at each location
advect_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to mean advective speed
starting_variance = NaN(length(ostrme),length(ostrme(1).lon)); %variance of salinity at beginning of streamline

maxlag = 96; %months
for ii = streamlines' %1:length(ostrme) %run through all streamlines on this surface
    disp(ii)
    if length(ostrme(ii).lat) ~= length(ostrme(1).lat)
        continue
    end
    latpoint(ii,:) = ostrme(ii).lat;
    lonpoint(ii,:) = ostrme(ii).lon;
    distkm(ii,:) = ostrme(ii).cdxkm;
    
    lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %longitude values of streamline
    lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat)); %latitude values of streamline
    if length(lat_nums) - length(lon_nums) == 1 %sometimes length of these are 1 off?
        lat_nums = lat_nums(1:(end-1));
    elseif length(lon_nums) - length(lat_nums) == 1
        lon_nums = lon_nums(1:(end-1));
    end
    cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %cumulative distance in km

    mean_pr(ii,1:length(lon_nums)) = interp2(a_xi,a_yi,sig_a_mean_pr',lon_nums,lat_nums,'linear');
     
    if length(lon_nums) < 21
        continue
    end
    
    %use mean advective speed 
    days_nums = ostrme(ii).days(~isnan(ostrme(ii).days)); %ecco streamline days since start, based on mean advective speed
    %for some reason ecco days not always same length as other variables...
    if length(days_nums) < length(cdxkm_nums)
        %lastdays = days_nums(end)-days_nums(end-1); %add one more day entry with same spacing as previous time step
        days_nums = [days_nums;zeros((length(cdxkm_nums)-length(days_nums)),1)];
    end
    
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

    %find the first point along the ecco streamline that has no NaNs, then go a
    %little past that
    sp = 1;
    %get an error when sp longer than 1st dim of e_strmln_spice_rm_seasonal
    stop = size(a_strmln_spice,1) - 20;
    while any(isnan(a_strmln_spice(sp,:))) && sp < stop
        sp = sp + 1;
    end
    outcrop = sp;
    if sp == stop
        continue
    end

    startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
    starting_spice = a_strmln_spice(startingpoint,:); %or strmln_spice_an
    starting_variance(ii,:) = var(starting_spice,'omitnan'); %calculate variance at beginning of streamline 
    a_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    a_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
    for pp = 1:length(lon_nums) %for each point along the streamline
        downstream_spice = a_strmln_spice(pp,:); %or strmln_spice_an (anomaly from overall mean)
        [r,lags] = xcorr(starting_spice,downstream_spice,maxlag,'coeff');
        a_R(:,pp) = r;
        [max_corr(ii,pp),ia] = max(r(lags <= 0)); %% make sure lag is negative
        max_lag(ii,pp) = lags(ia);
        if sum(~isnan(r)) == 0
            max_corr(ii,pp) = NaN;
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

    a_midx = NaN(1,length(lon_nums));
    for pp = 1:length(lon_nums)
        diff = NaN(length(lags),1);
        for ll = 1:length(lags)
            diff(ll) = abs(lags(ll) - months_nums(pp));
        end
        [~,a_midx(pp)] = min(diff);
        advect_corr(ii,pp) = a_R(length(lags)-a_midx(pp)+1,pp);
        advect_lag(ii,pp) = lags(a_midx(pp)); 
    end
end

%turn (streamline,point) arrays into 1D vectors
latpoint = reshape(latpoint,[length(ostrme)*3000,1]);
lonpoint = reshape(lonpoint,[length(ostrme)*3000,1]);
distkm = reshape(distkm,[length(ostrme)*3000,1]);
mean_pr = reshape(mean_pr,[length(ostrme)*3000,1]);
max_lag =  reshape(max_lag,[length(ostrme)*3000,1]);
advect_lag = reshape(advect_lag,[length(ostrme)*3000,1]);
max_corr = reshape(max_corr,[length(ostrme)*3000,1]);
advect_corr = reshape(advect_corr,[length(ostrme)*3000,1]);
starting_variance = reshape(starting_variance,[length(ostrme)*3000,1]);

%remove points with no data!
keep = ~isnan(latpoint) & ~isnan(lonpoint);
latpoint = latpoint(keep);
lonpoint = lonpoint(keep);
distkm = distkm(keep);
mean_pr = mean_pr(keep);
max_lag =  max_lag(keep);
advect_lag = advect_lag(keep);
max_corr = max_corr(keep);
advect_corr = advect_corr(keep);
starting_variance = starting_variance(keep);

%save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/argo_down_streamline_persistence_sig',sigma1,'.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');
save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/argo_down_streamline_persistence_sig',sigma1,'_tunnels.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');

%------------------------------------------------------------------------------------------------------
% Argo: for each streamline, calculate persistence of PV anomalies at each point

latpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of latitude coordinates
lonpoint = NaN(length(ostrme),length(ostrme(1).lon)); %list of longitude coordinates
distkm = NaN(length(ostrme),length(ostrme(1).lon)); %list of distances in km from starting point at surface outcrop
mean_pr = NaN(length(ostrme),length(ostrme(1).lon)); %list of mean pressure values at each point
max_corr = NaN(length(ostrme),length(ostrme(1).lon)); %max correlation at each location along streamline
advect_corr = NaN(length(ostrme),length(ostrme(1).lon)); %correlation following mean advective speed
max_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to max correlation at each location
advect_lag = NaN(length(ostrme),length(ostrme(1).lon)); %lag corresponding to mean advective speed
starting_variance = NaN(length(ostrme),length(ostrme(1).lon)); %variance of salinity at beginning of streamline

rho_ref = 1000;

maxlag = 96; %months

%get Argo data on this surface
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

for ii = streamlines' %1:length(ostrme) %run through all streamlines on this surface
    disp(ii)
    if length(ostrme(ii).lat) ~= length(ostrme(1).lat)
        continue
    end
    latpoint(ii,:) = ostrme(ii).lat;
    lonpoint(ii,:) = ostrme(ii).lon;
    distkm(ii,:) = ostrme(ii).cdxkm;
    
    lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %longitude values of streamline
    lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat)); %latitude values of streamline
    if length(lat_nums) - length(lon_nums) == 1 %sometimes length of these are 1 off?
        lat_nums = lat_nums(1:(end-1));
    elseif length(lon_nums) - length(lat_nums) == 1
        lon_nums = lon_nums(1:(end-1));
    end
    cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %cumulative distance in km
    mean_pr(ii,1:length(lon_nums)) = interp2(a_xi,a_yi,sig_a_mean_pr',lon_nums,lat_nums,'linear');
    
    if length(lon_nums) < 21
       continue
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
    
    a_strmln_pv = NaN(length(lon_nums),length(a_yrgrid_trunc)); %pv anomaly along streamline

    for it = 1:length(a_yrgrid_trunc)
        % fill out array of pv anomaly from mean, along ecco streamline
        a_strmln_pv(:,it) = interp2(a_xi,a_yi,squeeze(a_pv(:,:,it+12))',lon_nums,lat_nums,'linear'); %look into to see if somehow remaining nans where there shouldn't be
    end

    a_strmln_pv = a_strmln_pv(:,2:end); %for some reason the first month is NaNs

    % ARGO: find lagged correlation coefficients from point near outcrop (past
    % seasonal cycle) and points farther down streamline
    
    %find the first point along the ecco streamline that has no NaNs, then go a
    %little past that
    sp = 1;
    %get an error when sp longer than 1st dim of e_strmln_spice_rm_seasonal
    stop = size(a_strmln_pv,1) - 20;
    while any(isnan(a_strmln_pv(sp,:))) && sp < stop
        sp = sp + 1;
    end
    if sp == stop
        continue
    end

    startingpoint = sp+20; %how far down the streamline do you want to start? not km, just index of point
    a_starting_pv = a_strmln_pv(startingpoint,:); %or strmln_spice_an
    starting_variance(ii,:) = var(a_starting_pv,'omitnan'); %calculate variance at beginning of streamline 
    a_R = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    a_Lags = NaN(2*maxlag+1,length(lon_nums)); %-startingpoint
    a_advect = NaN(1,length(lon_nums)); %save lags that have max correlation for each location
    a_corrpeak = NaN(1,length(lon_nums)); %value of correlation peak
    for pp = 1:length(lon_nums) %for each point along the streamline
        a_downstream_pv = a_strmln_pv(pp,:); %or strmln_spice_an (anomaly from overall mean)
        [r,lags] = xcorr(a_starting_pv,a_downstream_pv,maxlag,'coeff');
        a_R(:,pp) = r;
        a_Lags(:,pp) = lags;
        [max_corr(ii,pp),ia] =  max(r(lags <= 0)); %% make sure lag is negative
        max_lag(ii,pp) = lags(ia);
        if sum(~isnan(r)) == 0
            max_lag(ii,pp) = NaN;
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

    a_midx = NaN(1,length(lon_nums));
    for pp = 1:length(lon_nums)
        diff = NaN(length(lags),1);
        for ll = 1:length(lags)
            diff(ll) = abs(lags(ll) - months_nums(pp));
        end
        [~,a_midx(pp)] = min(diff);
        advect_corr(ii,pp) = a_R(length(lags)-a_midx(pp)+1,pp);
        advect_lag(ii,pp) = lags(a_midx(pp)); 
    end   

end

%turn (streamline,point) arrays into 1D vectors
latpoint = reshape(latpoint,[length(ostrme)*3000,1]);
lonpoint = reshape(lonpoint,[length(ostrme)*3000,1]);
distkm = reshape(distkm,[length(ostrme)*3000,1]);
mean_pr = reshape(mean_pr,[length(ostrme)*3000,1]);
max_lag =  reshape(max_lag,[length(ostrme)*3000,1]);
advect_lag = reshape(advect_lag,[length(ostrme)*3000,1]);
max_corr = reshape(max_corr,[length(ostrme)*3000,1]);
advect_corr = reshape(advect_corr,[length(ostrme)*3000,1]);
starting_variance = reshape(starting_variance,[length(ostrme)*3000,1]);

%remove points with no data!
keep = ~isnan(latpoint) & ~isnan(lonpoint);
latpoint = latpoint(keep);
lonpoint = lonpoint(keep);
distkm = distkm(keep);
mean_pr = mean_pr(keep);
max_lag =  max_lag(keep);
advect_lag = advect_lag(keep);
max_corr = max_corr(keep);
advect_corr = advect_corr(keep);
starting_variance = starting_variance(keep);

%save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/PV_argo_down_streamline_persistence_sig',sigma1,'.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');
save(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/PV_argo_down_streamline_persistence_sig',sigma1,'_tunnels.mat'],'sigma1','latpoint','lonpoint','distkm','mean_pr','max_corr','advect_corr','max_lag','advect_lag','starting_variance');


%% plot 4 panels: argo spice, ecco spice, argo PV, ecco PV
sigma1 = '30.60';

figure('Position',[10 10 1400 500])
lvl = 0:0.1:1;

xrange = [0 360]; %atlantic
%xrange = [100 290]; %pacific
yrange = [-40 40];

custommap = [212 183 247;137 178 245;59 217 101;247 132 17;219 33 9]./252;
%argo spice
subaxis(2,2,1,'SpacingVert',0.035,'SpacingHoriz',0.02)
load(['/home/chersh/SpiceAnomalies/persistence_maps/argo_down_streamline_persistence_sig',sigma1,'.mat']);
lonpoint(lonpoint < 0) = lonpoint(lonpoint < 0) + 360;
xlabel('longitude')
ylabel('latitude')
scatter(lonpoint,latpoint,6,advect_corr,'filled')
hold on
%cmocean('amp',5);
%cmap = clmap(co);
colormap(custommap)
xlim(xrange)
ylim(yrange)
caxis([0 1])
set(gca,'XTick',[]);
vax = axis;
text(vax(1)+40,vax(4)-6.5,'Argo/spice','FontSize',8)
box on
text(vax(1)+280,vax(4)+12,['\sigma_1 = ',sigma1,' kg m^{-3}'],'FontSize',11)
yticks([-40 -20 0 20 40])
yticklabels({'40S','20S','0','20N','40N'})
gebco('k')
    
%ecco spice
load(['/home/chersh/SpiceAnomalies/persistence_maps/ecco_down_streamline_persistence_sig',sigma1,'.mat']);
subaxis(2,2,2)
xlabel('longitude')
ylabel('latitude')
lonpoint(lonpoint < 0) = lonpoint(lonpoint < 0) + 360;
scatter(lonpoint,latpoint,6,advect_corr,'filled')
hold on
colormap(custommap)
%cmocean('amp',5);
%cmap = clmap(co);
xlim(xrange)
ylim(yrange)
caxis([0 1])
set(gca,'XTick',[],'YTick',[]);
vax = axis;
text(vax(1)+40,vax(4)-6.5,'ECCO/spice','FontSize',8)
gebco('k')
box on

% argo PV
load(['/home/chersh/SpiceAnomalies/persistence_maps/PV_argo_down_streamline_persistence_sig',sigma1,'.mat']);
subaxis(2,2,3)
xlabel('longitude')
ylabel('latitude')
lonpoint(lonpoint < 0) = lonpoint(lonpoint < 0) + 360;
scatter(lonpoint,latpoint,6,advect_corr,'filled')
hold on
%cmocean('amp',5);
%cmap = clmap(co);
colormap(custommap)
xlim(xrange)
ylim(yrange)
caxis([0 1])
vax = axis;
text(vax(1)+40,vax(4)-6.5,'Argo/PV','FontSize',8)
yticks([-40 -20 0 20 40])
yticklabels({'40S','20S','0','20N','40N'})
xticklabels({'0','50E','100E','150E','160W','110W','60W','10W'})
gebco('k')
box on

% ecco PV
load(['/home/chersh/SpiceAnomalies/persistence_maps/PV_ecco_down_streamline_persistence_sig',sigma1,'.mat']);
subaxis(2,2,4)
xlabel('longitude')
ylabel('latitude')
lonpoint(lonpoint < 0) = lonpoint(lonpoint < 0) + 360;
scatter(lonpoint,latpoint,6,advect_corr,'filled')
hold on
%cmocean('amp',5);
%cmap = clmap(co);
cmap = colormap(custommap);
xlim(xrange)
ylim(yrange)
caxis([0 1])
set(gca,'YTick',[]);
vax = axis;
text(vax(1)+40,vax(4)-6.5,'ECCO/PV','FontSize',8)
c1 = colorbar;
c1.Label.String = 'correlation';
c1.Orientation = 'Vertical';
c1.Position = [.92 .1 .02 .8];
xticklabels({'0','50E','100E','150E','160W','110W','60W','10W'})
gebco('k')
box on

%saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/argo_ecco_spice_PV_persistence_sig',sigma1,'.png'])
%saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/argo_ecco_spice_PV_persistence_sig',sigma1,'.fig'])

%% plot 4 panels of lag on different surfaces

sigma1levels = {'28.05','29.95','29.55','30.60'};
figure('Position',[10 10 1400 500])
lvl = 0:0.1:1;

xrange = [0 360]; %atlantic
%xrange = [100 290]; %pacific
yrange = [-40 40];

%lag
for ss = 1:length(sigma1levels)
    load(['/atalanta/home/chersh/SpiceAnomalies/persistence_maps/argo_down_streamline_persistence_sig',sigma1levels{ss},'.mat']);
    subaxis(2,2,ss,'SpacingVert',0.03,'SH',0.03)
    xlabel('longitude')
    ylabel('latitude')
    scatter(lonpoint,latpoint,6,advect_lag./12,'filled')
    hold on
    gebco('k')
    xlim(xrange)
    ylim(yrange)
    if ss == 4
        c1 = colorbar;
        c1.Label.String = 'years';
        c1.Orientation = 'Horizontal';
        c1.Position = [.3 .91 .4 .01];
    end
    if ss == 1
        set(gca,'XTick',[]);
    elseif ss == 2
        set(gca,'XTick',[],'YTick',[])
    elseif ss == 4
        set(gca,'YTick',[])
    end
    caxis([0 8])
    colormap(turbo(8))
    %cmocean('thermal',6);
    %cmap = clmap(18);
    vax = axis;
    text(vax(1)+40,vax(4)-6,['\sigma_1 = ',sigma1levels{ss}],'FontSize',8)
end

%sgtitle('advective time scale')

saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/ecco_advection_lag.png')
saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/ecco_advection_lag.fig')
%% convert to netcdf
out_nc = 'PV_ecco_down_streamline_persistence_sig30.60.nc';
nccreate(out_nc,'advect_corr','Dimensions',{'npoints',length(advect_corr)});
ncwrite(out_nc,'advect_corr',advect_corr);
nccreate(out_nc,'advect_lag','Dimensions',{'npoints',length(advect_corr)});
ncwrite(out_nc,'advect_lag',advect_lag);
nccreate(out_nc,'latpoint','Dimensions',{'npoints',length(advect_corr)});
ncwrite(out_nc,'latpoint',latpoint);
nccreate(out_nc,'lonpoint','Dimensions',{'npoints',length(advect_corr)});
ncwrite(out_nc,'lonpoint',lonpoint);


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



