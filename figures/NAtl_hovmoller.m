%% from compare_argo_ecco4r4_2023.m
%C Hersh, 2023

% plots spiciness and PV in Argo analysis and ECCO

%if on taris
%{
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/SpiceAnomalies
addpath /atalanta/home/chersh/SpiceAnomalies/strmat_sla_200rmax
addpath /atalanta/home/chersh/toolbox/cmocean
cd /atalanta/home/swijffels/work/argo/gridNSF
%}

addpath /home/swijffels/toolbox/seawater
addpath /home/swijffels/toolbox/csirolib
addpath /home/swijffels/work/argo/matlab
addpath /home/swijffels/toolbox/subaxis
addpath /home/chersh/SpiceAnomalies
addpath /home/chersh/SpiceAnomalies/strmat_sla_200rmax
addpath /home/chersh/toolbox/cmocean
cd /home/swijffels/work/argo/gridNSF

%% load in climatology 
if ~exist('clim')
    clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end 

aclim_sig1grid = clim.sig1grid;

a_mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,1:85)) + squeeze(clim.pr.c(1:4:end,21:4:621,1:85,1));

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end

%% WHICH SIGMA1 LEVEL DO YOU WANT TO LOOK AT?
%set(gcf,'PaperPositionMode','auto') -- saves figure as what you see on
%screen
sig1 = '30.60';
[~,aclim_isig]=min(abs(aclim_sig1grid - str2double(sig1)));

%load in the streamline mat file for a particular density level
load(['ECCO4r4_strm_ventilated_',sig1,'.mat']);

%WHICH STREAMLINE DO YOU WANT TO LOOK AT?
ii = 142; %choose by visually identifying ventilation regions on streamline plots

% get ecco data
e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');
e_sig1grid = e.sig1grid; 
[~,e_isig]=min(abs(e_sig1grid - str2double(sig1))); %convert sigma1 value to index of sigma grid

%  make Pacific centric
if min(e.xi) < 0
    xo = e.xi;
    xi = rem(xo+360,360);
    
    [ecco.xi,ixs]=sort(xi);
    
    junk = squeeze(e.si(:,:,e_isig,:));
    ecco.si = junk(ixs,:,:);
    junk = squeeze(e.pri(:,:,e_isig,:));
    ecco.pri = junk(ixs,:,:);
end

% load in Argo data analysis
if ~exist('sa_fwa')
    %a = matfile('/atalanta/home/chersh/SpiceAnomalies/despike_argo_data.mat');
    a = matfile('/home/swijffels/work/argo/gridNSF/gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
    a_sig1grid = a.sig1grid;
    [~,a_isig]=min(abs(a_sig1grid - str2double(sig1))); %convert sigma1 value to index of sigma grid

    sa = squeeze(a.sa_fwa(:,:,a_isig,:)); %fwa is 3-month moving window
    pr = squeeze(a.pr_fwa(:,:,a_isig,:));
    dpr = squeeze(a.dpr_fwa(:,:,a_isig,:));
    dpr_below = squeeze(a.dpr_fwa(:,:,a_isig+1,:)); %layer thickness between this level and one below
end

% trim sigma analysis using sfc values -- question: should do this for ecco as well??
%for it = 1:nt 
%    z0 = squeeze(sfc.sig1(:,:,1,it));
%    for is = isig %find(sig1grid<= nanmax(z0(:)))' 
%        imask = NaN(size(z0));
%        imask(z0 <= sig1grid(is)) = 1;
%        sa(:,:,it) = squeeze(sa(:,:,it)).*imask;
%        pr(:,:,it) = squeeze(pr(:,:,it)).*imask;
%        dpr(:,:,it) = squeeze(dpr(:,:,it)).*imask;
%    end
%end

%% load in and process all data to plot

%SPICE
%get mean argo salinity climatology
mean_sa = squeeze(clim.sa.m(:,:,aclim_isig)) + squeeze(clim.sa.c(:,:,aclim_isig,1));
mean_sa = mean_sa(1:4:end,21:4:621);
mean_sa(squeeze(a_mean_pr(:,:,a_isig)) <= clim_mldminmax.mldmax(1:4:end,21:4:621)) = NaN;
mean_sa(str2double(sig1) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621)) = NaN;

nstrm = length(ostrme); %how many streamlines are available for this sigma1 level?

yi = a.yi;
xi = a.xi;
yrgrid = a.yrgrid;
[nx,ny,nt]=size(sa);

%pick some distances to mark on map
markdist = 0:2000:6000;

cdxkm_nums = ostrme(ii).cdxkm(~isnan(ostrme(ii).cdxkm)); %some streamlines have NaNs at the end of cdxkm vector
days_nums = ostrme(ii).days(~isnan(ostrme(ii).days)); %get rid of NaNs at end of days vector

%for some reason ecco days not always same length as other variables...
if length(days_nums) < length(cdxkm_nums)
    %lastdays = dayse_nums(end)-dayse_nums(end-1); %add one more day entry with same spacing as previous time step
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

rossbyspeed = ostrme(ii).rossbyspeed(1:length(cdxkm_nums)); %add background U to rossby speed
rossby_days = ((dxkm*1000)./rossbyspeed)/(60*60*24); %calculate days from Rossby phase speed -- ecco path
rossby_days = cumsum(rossby_days);
advect_linex = cdxkm_nums;
advect_liney = 2004 + days_nums./365;
%speed of rossby wave propagation
rossby_linex = cdxkm_nums;
rossby_liney = 2004 + rossby_days./365;

lon_nums = ostrme(ii).lon(~isnan(ostrme(ii).lon)); %ecco
lat_nums = ostrme(ii).lat(~isnan(ostrme(ii).lat));
%found case where length(lon_nums) ~= length(lat_nums)
if length(lon_nums) > length(lat_nums)
    lon_nums = lon_nums(1:length(lat_nums)); %cut off extra point(s)
elseif length(lat_nums) > length(lon_nums)
    lat_nums = lat_nums(1:length(lon_nums));
end

markdistidx = NaN(1,length(markdist));
for dd = 1:length(markdist)
    [~,markdistidx(dd)] = min(abs(markdist(dd)-cdxkm_nums));
end

res = 1; %what resolution do you actually want to plot?
distpath = cdxkm_nums(1:res:end);
xpath = lon_nums(1:res:end);
ypath = lat_nums(1:res:end);

%argo on ecco streamline
zp = deal(NaN(length(xpath),nt));
% pull out data along this path - Argo
for it = 1:nt
    %argo on ecco streamline
    zp(:,it) = interp2(xi,yi,squeeze(sa(:,:,it))',xpath,ypath,'linear');
end

argo_spice_data = zp;

% pull out data along the streamline - ecco
zpe = NaN(length(xpath),length(e.yrgrid));
mean_sa_e = nanmean(squeeze(e.si(:,:,e_isig,:)),3); %pull out mean ecco salinity on this sigma1 level
zpme  = interp2(e.xi,e.yi,squeeze(mean_sa_e)',xpath,ypath,'linear');

for it = 1:length(e.yrgrid)
   %anomaly from ecco mean, along ecco streamline
   zpe(:,it) = interp2(e.xi,e.yi,squeeze(e.si(:,:,e_isig,it))',xpath,ypath,'linear')-zpme;
end

for ll = 1:size(zpe,1) %remove seasonal cycle to match argo
    if all(isnan(squeeze(zpe(ll,:))))
        continue
    end
    [seasonal_anom,seasonal_mean,seasonal_coef,seasonal_m] = seasonal(e.yrgrid,squeeze(zpe(ll,:)));
    zpe(ll,:) = zpe(ll,:) - seasonal_mean';
end

ecco_spice_data = zpe;

%PV
%-------------------------------------------------------------------------------------------
rho_ref = 1000; %reference density for calculating pv

%layer thickness climatology (pressure diff between sigma layer and the one above it)
dpclim2 = squeeze(clim.dpr.m(:,:,aclim_isig)) + squeeze(clim.dpr.c(:,:,aclim_isig,1));
%layer thickness climatology (pressure diff between sigma layer and the one below it)
dpclim1 = squeeze(clim.dpr.m(:,:,aclim_isig+1)) + squeeze(clim.dpr.c(:,:,aclim_isig+1,1));

fcor = calc_fcor(clim.xi,clim.yi); %coriolis parameter
drho2 = a_sig1grid(a_isig) - a_sig1grid(a_isig-1); %change in density between this layer and the one above it
drho1 = a_sig1grid(a_isig+1) - a_sig1grid(a_isig); %change in density between this layer and the one below it

drhodpclim = (drho2./dpclim2 + drho1./dpclim1)/2;
pvclim = (1/rho_ref)*fcor.*drhodpclim; %pv climatological mean overall
pv_anom =  NaN(size(dpr)); %calculate pv anomaly for argo
for it=1:nt
    dpr2 = dpr(:,:,it); %anomaly of layer thickness between this layer and one above it
    dpr1 = dpr_below(:,:,it); %anomaly of layer thickness between this layer and one below it
    drhodp_tot = (drho2./(dpclim2(1:4:end,21:4:621) + dpr2) + drho1./(dpclim1(1:4:end,21:4:621) + dpr1))/2;
    pv_anom(:,:,it) = (1/rho_ref)*fcor(1:4:end,21:4:621).*drhodp_tot - pvclim(1:4:end,21:4:621);
end

pvclim = pvclim(1:4:end,21:4:621);
pvclim(squeeze(a_mean_pr(:,:,a_isig)) <= clim_mldminmax.mldmax(1:4:end,21:4:621)) = NaN;
pvclim(str2double(sig1) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621)) = NaN;

%argo on ecco streamline
zp = NaN(length(xpath),nt); %initialize matrix for argo anomalies along streamline

% pull out data along this path - Argo (but use ECCO streamline!)
for it = 1:nt
    %argo on ecco streamline
    zp(:,it) = interp2(xi,yi,squeeze(pv_anom(:,:,it))',xpath,ypath,'linear');
end

argo_pv_data = zp;

% ECCO: pull out data along the streamline
% use the ecco streamline
zpe = NaN(length(xpath),length(e.yrgrid));
drhodr_ecco = squeeze(e.drhodri(:,:,e_isig,:));
drhodr_ecco_mean = nanmean(drhodr_ecco,3);

%recalculate fcor with new dimensions
fcor = calc_fcor(e.xi,e.yi); %coriolis parameter
yi_e = e.yi;
xi_e = e.xi;

pvclim_ecco = (1/rho_ref)*fcor.*drhodr_ecco_mean;
zpme = interp2(xi_e,yi_e,squeeze(pvclim_ecco)',xpath,ypath,'linear');

for it = 1:length(e.yrgrid)
   % anomaly from ecco mean, along ecco streamline
   zpe(:,it) = interp2(xi_e,yi_e,(((1/rho_ref)*fcor.*drhodr_ecco(:,:,it)))',xpath,ypath,'linear')-zpme;
end
zpe = -zpe; %because ECCO DRHODR variable is negative (r points up)

for ll = 1:size(zpe,1) %remove seasonal cycle to match argo
    if sum(isnan(squeeze(zpe(ll,:)))) > 295
        continue
    end
    [seasonal_anom,seasonal_mean,seasonal_coef,seasonal_m] = seasonal(e.yrgrid,squeeze(zpe(ll,:)));
    zpe(ll,:) = zpe(ll,:) - seasonal_mean';
end
ecco_pv_data = zpe;

%% mean maps and Hovmoller plots

figure('Position',[10 -100 1600 800])

% plot mean argo salinity with streamlines overlaid
[sprc]=prctile(mean_sa(:)',[5,99]);
ax1 = subaxis(1,3,1,'SH',0.02);
cl = [sprc(1)+ [0:0.025:1]*(sprc(2) - sprc(1))];
contourf(xi,yi,mean_sa',[20,cl,40],'linecolor','flat') 
hold on
caxis(sprc);
gebco('k')
cmap = cmocean('haline');
colormap(ax1,cmap)
c = colorbar('location','southoutside');
c.Label.String = 'psu';
c.Label.FontSize = 10;
c.FontSize = 10;
title(['{\sigma}_1 = ',num2str(ostrme(ii).sigma1,4)],'FontSize',16);
axis([260,350,-40,40])
yticklabels({'40S','30S','20S','10S','0','10N','20N','30N','40N'})
xticks([260 280 300 320 340])
xticklabels({'100W','80W','60W','40W','20W'})
ax1.XAxis.FontSize = 8;
ax1.YAxis.FontSize = 8;

%plot streamlines on map
hold on
xpath_plot = xpath(distpath <= 6000);
ypath_plot = ypath(distpath <= 6000);
for dd = 1:length(markdist)
    plot(xpath(markdistidx(dd)),ypath(markdistidx(dd)),'r*','LineWidth',2)
end
plot(xpath(1),ypath(1),'r*','LineWidth',3) %streamline start
plot(xpath_plot,ypath_plot,'r-','LineWidth',2) %ecco streamline

%plot argo hovmoller - clmap(24)?
ax2 = subaxis(1,3,2);
%argo on ecco streamline
[sprc]=prctile(argo_spice_data(:)',[5,95]);
pcolor(distpath,yrgrid,argo_spice_data'),shading flat
hold on
plot(rossby_linex,2004+rem(rossby_liney-2004,15.5),'k:','Linewidth',4);
%advect_linex(328:end) = NaN;
%advect_liney(328:end) = NaN;
plot(advect_linex,2004+rem(advect_liney-2004,15.5),'k','Linewidth',4);
xlim([0 6000])
cmap = cmocean('balance');
colormap(ax2,cmap)
caxis([-1,1]*max(abs(sprc(:))))
va = axis;
axis([va(1:2),2000,2020])
hold on
plot(va(1:2),[2004,2004],'k-')
plot(va(1:2),[2016,2016],'k-')
c3 = colorbar('FontSize',10,'Location','southoutside');
c3.Label.String = 'psu';
%c3.Position = [0.068 0.1 .02 .375];
%c3.Orientation = 'Vertical';
xlabel('km','FontSize',12)
set(gca,'YTick',[])
title('Argo','FontSize',14);

%plot ecco hovmoller
ax3 = subaxis(1,3,3);
%if using ecco streamline
pcolor(distpath,e.yrgrid,ecco_spice_data'),shading flat
hold on
plot(rossby_linex,2004+rem(rossby_liney-2004,15.5),'k:','Linewidth',4); %days not filled out
plot(advect_linex,2004+rem(advect_liney-2004,15.5),'k','Linewidth',4);
caxis([-1,1]*max(abs(sprc(:))))
cmap = cmocean('balance');
colormap(ax3,cmap)
xlim([0 6000])
va = axis;
axis([va(1:2),2000,2020])
ylabel('Year','FontSize',12)
xlabel('km','FontSize',12)
set(gca,'YAxisLocation','right')
title('ECCO','FontSize',14);
hold on
plot(va(1:2),[2004,2004],'k-')
plot(va(1:2),[2016,2016],'k-')

%saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_hovmoller_spice.png')
%saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_hovmoller_spice.fig')

% PV DRHODR -mean  maps and pull out float pathways

%plot mean argo PV and use to make pathways

figure('Position',[10 -100 1600 800])

ax1 = subaxis(1,3,1,'SH',0.02);
[sprc]=prctile(abs(pvclim(:)'),[5,95]);
cl = [sprc(1)+ [0:0.025:1]*(sprc(2) - sprc(1))];
contourf(xi,yi,abs(pvclim)',[-1,cl,1],'linecolor','flat') 
hold on
caxis([sprc(1) sprc(2)]);
gebco('k')
c1 = colorbar;
cmocean('haline');
c1.Label.String = 'kg s^{-1} m^{-4}';
c1.Label.FontSize = 10;
c1.Location = 'southoutside';
title(['{\sigma}_1 = ',num2str(ostrme(ii).sigma1,4)],'FontSize',16);
axis([260,350,-40,40])
yticklabels({'40S','30S','20S','10S','0','10N','20N','30N','40N'})
xticks([260 280 300 320 340])
xticklabels({'100W','80W','60W','40W','20W'})
ax1.XAxis.FontSize = 8;
ax1.YAxis.FontSize = 8;
%plot streamlines on map
hold on
for dd = 1:length(markdist)
    plot(xpath(markdistidx(dd)),ypath(markdistidx(dd)),'r*','LineWidth',2)
end
plot(xpath(1),ypath(1),'r*','LineWidth',3) %streamline start
plot(xpath_plot,ypath_plot,'r-','LineWidth',2) %ecco streamline

[sprc]=prctile(argo_pv_data(:)',[5,95]); %get bounds for colorbar

%plot argo hovmoller
ax5 = subaxis(1,3,2);
%argo on ecco streamline
pcolor(distpath,yrgrid,argo_pv_data'),shading flat
hold on
plot(rossby_linex,2004+rem(rossby_liney-2004,15.5),'k:','Linewidth',4);
plot(advect_linex,2004+rem(advect_liney-2004,15.5),'k','Linewidth',4);
caxis([-1,1]*max(abs(sprc(:))))
colormap(ax5,cmap)
va = axis;
axis([va(1:2),2000,2020])
xlim([0 6000])
hold on
plot(va(1:2),[2004,2004],'k-')
plot(va(1:2),[2016,2016],'k-')
c4 = colorbar('FontSize',10,'Location','southoutside');
c4.Label.String = 'kg s^{-1} m^{-4}';
%c4.Position = [0.068 0.1 .02 .375];
%c4.Orientation = 'Vertical';
xlabel('km','FontSize',12)
set(gca,'YTick',[])
title('Argo','FontSize',14);

%plot ecco hovmoller
ax6 = subaxis(1,3,3);
%using ecco streamline
pcolor(distpath,e.yrgrid,ecco_pv_data'),shading flat
hold on
plot(rossby_linex,2004+rem(rossby_liney-2004,15.5),'k:','Linewidth',4);
plot(advect_linex,2004+rem(advect_liney-2004,15.5),'k','Linewidth',4);
caxis([-1,1]*max(abs(sprc(:))))
xlim([0 6000])
va = axis;
axis([va(1:2),2000,2020])
colormap(ax6,cmap)
ylabel('Year','FontSize',12)
xlabel('km','FontSize',12)
set(gca,'YAxisLocation','right')
title('ECCO','FontSize',14);
hold on
plot(va(1:2),[2004,2004],'k-')
plot(va(1:2),[2016,2016],'k-')

%saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_hovmoller_PV.png')
%saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/NAtl_hovmoller_PV.fig')

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