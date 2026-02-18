%% add paths
% look at a horizontal/depth snapshot of spice anomalies 
%if on taris
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/SpiceAnomalies
addpath /atalanta/home/chersh/SpiceAnomalies/strmat_sla_200rmax
addpath /atalanta/home/chersh/toolbox/cmocean
cd /atalanta/home/swijffels/work/argo/gridNSF

%% load in climatology 
if ~exist('clim')
    clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end 

aclim_sig1grid = clim.sig1grid;

a_mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,1:85)) + squeeze(clim.pr.c(1:4:end,21:4:621,1:85,1));

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/atalanta/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end

%% choose streamline and month pull out data
sig1 = '28.95';
[~,aclim_isig]=min(abs(aclim_sig1grid - str2double(sig1)));


nmonth = 100; %argo starts in 2003
e_nmonth = [nmonth+119,nmonth+120]; %ecco starts in 1993 and is half a month offset from argo;
%average month before and after

%load in the streamline mat file for a particular density level
load(['ECCO4r4_strm_ventilated_',sig1,'.mat']);

%WHICH STREAMLINE DO YOU WANT TO LOOK AT?
ii = 52; %choose by visually identifying ventilation regions on streamline plots

% get ecco data
e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');
e_sig1grid = e.sig1grid; 
[~,e_isig]=min(abs(e_sig1grid - str2double(sig1))); %convert sigma1 value to index of sigma grid

%  make Pacific centric
if min(e.xi) < 0
    xo = e.xi;
    xi = rem(xo+360,360);
    
    [ecco.xi,ixs]=sort(xi);
    
    junk = squeeze(e.si(:,:,:,e_nmonth));
    ecco.si = junk(ixs,:,:);
    junk = squeeze(e.pri(:,:,:,e_nmonth));
    ecco.pri = junk(ixs,:,:);
end

% load in Argo data analysis
if ~exist('sa_fwa')
    %a = matfile('/atalanta/home/chersh/SpiceAnomalies/despike_argo_data.mat');
    a = matfile('/atalanta/home/swijffels/work/argo/gridNSF/gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
    a_sig1grid = a.sig1grid;
    [~,a_isig]=min(abs(a_sig1grid - str2double(sig1))); %convert sigma1 value to index of sigma grid

    sa = squeeze(a.sa_fwa(:,:,:,nmonth)); %fwa is 3-month moving window
    pr = squeeze(a.pr_fwa(:,:,:,nmonth));
end

%% load in and process all data to plot

%SPICE
%get mean argo salinity climatology
mean_sa = squeeze(clim.sa.m(:,:,:)) + squeeze(clim.sa.c(:,:,:,1));
mean_sa = mean_sa(1:4:end,21:4:621,:);
%mean_sa(squeeze(a_mean_pr(:,:,:)) <= clim_mldminmax.mldmax(1:4:end,21:4:621)) = NaN;
%mean_sa(str2double(sig1) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621)) = NaN;

yi = a.yi;
xi = a.xi;
[nx,ny,nsig]=size(sa);

cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %some streamlines have NaNs at the end of cdxkm vector
days_nums = ostrme(ii).days(~isnan(ostrme(ii).days)); %get rid of NaNs at end of days vector

%for some reason ecco days not always same length as other variables...
if length(days_nums) < length(cdxkm_nums)
    days_nums = [days_nums;zeros((length(cdxkm_nums)-length(days_nums)),1)];
end

if any(ostrme(ii).speed == inf)
    ostrme(ii).speed(ostrme(ii).speed==inf) = NaN; %set infinity values to NaN for now
end

dxkm = NaN(size(cdxkm_nums)); %find NON cumulative distances between streamline points -- ecco
for ll = 1:(length(cdxkm_nums)-1)
    dxkm(ll) = cdxkm_nums(ll+1) - cdxkm_nums(ll);
end
dxkm(end) = dxkm(end-1); %just set final distance to distance at penultimate point

lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %ecco
lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat));
%found case where length(lon_nums) ~= length(lat_nums)
if length(lon_nums) > length(lat_nums)
    lon_nums = lon_nums(1:length(lat_nums)); %cut off extra point(s)
elseif length(lat_nums) > length(lon_nums)
    lat_nums = lat_nums(1:length(lon_nums));
end

res = 1; %what resolution do you actually want to plot?
distpath = cdxkm_nums(1:res:end);
xpath = lon_nums(1:res:end);
ypath = lat_nums(1:res:end);

%argo on ecco streamline
zp = deal(NaN(length(xpath),nsig));
% pull out data along this path - Argo
for is = 1:nsig
    %argo on ecco streamline
    zp(:,is) = interp2(xi,yi,squeeze(sa(:,:,is))',xpath,ypath,'linear');
end

argo_spice_data = zp;

% pull out data along the streamline - ecco
zpe1 = NaN(length(xpath),length(e.sig1grid)); %ecco salinity anomaly at timestep before
zpe2 = NaN(length(xpath),length(e.sig1grid)); %ecco salinity anomaly at timestep after

for is = 1:length(e.sig1grid) %for each sigma1 level
   %anomaly from ecco mean, along ecco streamline
    mean_sa_e = nanmean(squeeze(e.si(:,:,is,:)),3); %pull out mean ecco salinity at this sigma1 level
    zpme  = interp2(e.xi,e.yi,squeeze(mean_sa_e)',xpath,ypath,'linear');

    zpe1(:,is) = interp2(e.xi,e.yi,squeeze(e.si(:,:,is,e_nmonth(1)))',xpath,ypath,'linear')-zpme;
    zpe2(:,is) = interp2(e.xi,e.yi,squeeze(e.si(:,:,is,e_nmonth(2)))',xpath,ypath,'linear')-zpme;
end

for ll = 1:size(zpe1,1) %remove seasonal cycle to match argo
    if all(isnan(squeeze(zpe1(ll,:))))
        continue
    end
    [seasonal_anom,seasonal_mean,seasonal_coef,seasonal_m] = seasonal(e.yrgrid,squeeze(zpe1(ll,:)));
    zpe1(ll,:) = zpe1(ll,:) - seasonal_mean';
end

for ll = 1:size(zpe2,1) %remove seasonal cycle to match argo
    if all(isnan(squeeze(zpe2(ll,:))))
        continue
    end
    [seasonal_anom,seasonal_mean,seasonal_coef,seasonal_m] = seasonal(e.yrgrid,squeeze(zpe2(ll,:)));
    zpe2(ll,:) = zpe2(ll,:) - seasonal_mean';
end

ecco_spice_data = (zpe1 + zpe2)/2;
