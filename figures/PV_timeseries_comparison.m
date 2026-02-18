addpath /atalanta/home/swijffels/toolbox/susan
addpath /atalanta/home/swijffels/toolbox/eez
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/toolbox
cd /atalanta/home/swijffels/work/argo/gridNSF

%change this to plot PV instead of salinity
%% load in argo data
if ~exist('sa_fwa')
anom = matfile('gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
end
yrgrid = anom.yrgrid;
a_xi = anom.xi;
a_yi = anom.yi;
sig1grid = anom.sig1grid;

%% get climatology:
if ~exist('clim')
%clim = load('gridonSig__climatology_12_2019.mat');
clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end

% make mean sigma depths:
zmean = squeeze(clim.pr.m) + squeeze(clim.pr.c(:,:,:,1));

%% get argo input data and R&G data
dat = load('Sigma1_anomalies_Argo_CORA_sla_superobs_1998_2020.mat','decyr','sa_sig','lat','lon','pa_sig','dpa_sig','sig1grid');

rg = matfile('/atalanta/home/chersh/SpiceAnomalies/RandG_on_sigma1.mat'); %Roemmich and Gilson gridded data, interpolated onto sigma1
rg_yrgrid = rg.yrgrid;
rg_xi = rg.xi;
rg_yi = rg.yi;
select_sig1grid = rg.select_sig1grid;
%% get ECCO data
e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');

e_xi = e.xi;
e_yi = e.yi;
e_yrgrid = e.yrgrid;
e_sig1grid = e.sig1grid;

%find overlap in time between Argo and ECCO
%Argo start is start bound, ECCO end is end bound
e_tidx = e_yrgrid > yrgrid(1);
a_tidx = yrgrid < e_yrgrid(end);

%% time series PV
addpath /atalanta/home/chersh/toolbox

xplot = [220,320,220,80];
yplot = [20,20,-20,-10];
is = [38,44,34,40]; %indices of sigma1 levels
name = {'North Pacific','North Atlantic','South Pacific','South Indian'};
positions = [0.08 0.57;0.53 0.7;0.08 0.13;0.57 0.13];

rho_ref = 1000;

figure('Position',[10 10 1600 900])
    
for ll = 1:length(xplot)

    ii = abs(dat.lon - xplot(ll))<2 & abs(dat.lat - yplot(ll))<1.5;

    disp(is(ll))
    listx = 1:length(a_xi); %this is stupid
    listy = 1:length(a_yi);
    idxx = listx(a_xi==xplot(ll));
    idxy = listy(a_yi==yplot(ll));
    
    % CALCULATE ARGO PV
    a_dpr = squeeze(anom.dpr_fwa(:,:,is(ll),:)); %layer thickness -- could just take out single timeseries here for argo (ecco needs to interpolate?)
    a_dpr_below = squeeze(anom.dpr_fwa(:,:,is(ll)+1,:)); %layer thickness between this level and one below
    
    %layer thickness climatology (pressure diff between sigma layer and the one above it)
    dpclim2 = squeeze(clim.dpr.m(:,:,is(ll))) + squeeze(clim.dpr.c(:,:,is(ll),1));
    %layer thickness climatology (pressure diff between sigma layer and the one below it)
    dpclim1 = squeeze(clim.dpr.m(:,:,is(ll)+1)) + squeeze(clim.dpr.c(:,:,is(ll)+1,1));
    
    %calculate fcor with clim dimensions
    fcor = calc_fcor(clim.xi,clim.yi);
    
    drho2 = sig1grid(is(ll)) - sig1grid(is(ll)-1); %change in density between this layer and the one above it
    drho1 = sig1grid(is(ll)+1) - sig1grid(is(ll)); %change in density between this layer and the one below it
    
    drhodpclim = (drho2./dpclim2 + drho1./dpclim1)/2;
    pvclim = (1/rho_ref)*fcor.*drhodpclim; %pv climatological mean overall
    a_pv =  NaN(length(a_xi),length(a_yi),length(yrgrid)); %calculate pv anomaly for argo

    for it=1:length(yrgrid)
        dpr2 = a_dpr(:,:,it); %anomaly of layer thickness between this layer and one above it
        dpr1 = a_dpr_below(:,:,it); %anomaly of layer thickness between this layer and one below it
        drhodp_tot = (drho2./(dpclim2(1:4:end,21:4:621) + dpr2) + drho1./(dpclim1(1:4:end,21:4:621) + dpr1))/2;
        a_pv(:,:,it) = (1/rho_ref)*fcor(1:4:end,21:4:621).*drhodp_tot - pvclim(1:4:end,21:4:621);
    end 
    a_pv = permute(a_pv,[2,1,3]);
    
    zz = squeeze(a_pv(idxx,idxy,is(ll),:));
    ze = squeeze(anom.dpr_fwe(idxx,idxy,is(ll),:));
    
    % CALCULATE ECCO PV
    fcor = calc_fcor(e_xi,e_yi);
    drhodr_sig = squeeze(e.drhodri(:,:,is(ll)-13,:)); %pull out drho/dr on this sigma surface
    e_mean_drhodr = nanmean(drhodr_sig,3);
    e_mean_pv = (1/rho_ref)*fcor.*e_mean_drhodr;
    e_pv = -1*((1/rho_ref)*fcor.*drhodr_sig - e_mean_pv); %times -1 b/c drhodri is negative
    
    [e_XI,e_YI] = meshgrid(e_xi,e_yi);
    [a_XI,a_YI] = meshgrid(a_xi,a_yi);
    
    e_pv_sparse = NaN(length(a_yi),length(a_xi),length(e_yrgrid));
    for tt = 1:length(e_yrgrid)
        e_pv_sparse(:,:,tt) = interp2(e_XI,e_YI,e_pv(:,:,tt)',a_XI,a_YI);
    end
    
    % get ecco timeseries
    [~,idxx] = min(abs(e_xi - xplot(ll)));
    [~,idxy] = min(abs(e_yi - yplot(ll)));
    e_tseries = squeeze(e_pv_sparse(idxx,idxy,is(ll)-13,:));
    e_tseries = e_tseries(e_tidx); %restrict to overlap with argo
    %now re-set the mean to zero now that I've restricted the
    %time series
    e_tseries = e_tseries - mean(e_tseries,'omitnan');
    e_seasonal_cycle = zeros(size(e_tseries));
    new_e_yrgrid = e_yrgrid(e_tidx); % restrict yrgrid to overlap with ecco
    
    %remove seasonal from ECCO  
    %calculate seasonal cycle for this timeseries
    if sum(~isnan(e_tseries)) < 5 %if no ecco data, don't bother making figure
        continue
    end
    [zanom,zmean,coef,zm] = seasonal(new_e_yrgrid,e_tseries);
    e_seasonal_cycle = zmean;
    e_timeseries = e_tseries - e_seasonal_cycle;
    
    % CALCULATE R&G PV
    [~,idxx] = min(abs(rg_xi - xplot(ll)));
    [~,idxy] = min(abs(rg_yi - yplot(ll)));
    is_rg = find(select_sig1grid == sig1grid(is(ll)));
    
    %calculate R&G layer thickness from sigma1 grid, pressure; then calculate pv
    %pres_mean_sig1, pres_anom_sig1
    %calculate R&G layer thickness climatology (pressure diff between sigma layer and the one above it)
    rg_dpr_clim = NaN(size(rg.pres_mean_sig1));
    for jj = 2:length(select_sig1grid)
        rg_dpr_clim(:,:,jj) =  rg.pres_mean_sig1(:,:,jj) - rg.pres_mean_sig1(:,:,jj-1);
    end
    rg_dpr_clim(:,:,1) = rg_dpr_clim(:,:,2); % just fill in the top sigma1 layer with layer below
    
    rg_dpr_anom = NaN(size(rg.pres_anom_sig1));
    for jj = 2:length(select_sig1grid) %for each sigma1 layer (36)
        for tt = 1:size(rg.pres_anom_sig1,4) %for each month
            rg_dpr_anom(:,:,jj,tt) = rg.pres_anom_sig1(:,:,jj,tt) - rg.pres_anom_sig1(:,:,jj-1,tt);
        end
    end
    rg_dpr_anom(:,:,1,:) = rg_dpr_anom(:,:,2,:); %just fill in the top sigma1 layer with layer below
    
    rg_dpr_clim_below = NaN(size(rg.pres_mean_sig1));
    for jj = 1:(length(select_sig1grid)-1)
        rg_dpr_clim_below(:,:,jj) = rg.pres_mean_sig1(:,:,jj+1) - rg.pres_mean_sig1(:,:,jj);
    end
    rg_dpr_clim_below(:,:,end) = rg_dpr_clim_below(:,:,2); % just fill in the bottom sigma1 layer with layer above
    
    rg_dpr_anom_below = NaN(size(rg.pres_anom_sig1));
    for jj = 1:(length(select_sig1grid)-1) %for each sigma1 layer (36)
        for tt = 1:size(rg.pres_anom_sig1,4) %for each month
            rg_dpr_anom_below(:,:,jj,tt) = rg.pres_anom_sig1(:,:,jj+1,tt) - rg.pres_anom_sig1(:,:,jj,tt);
        end
    end
    rg_dpr_anom_below(:,:,end,:) = rg_dpr_anom_below(:,:,end-1,:); %just fill in the bottom sigma1 layer with layer above
    
    fcor = calc_fcor(rg_xi,rg_yi);
    %drho2 and drho1 already calculated above for new argo product
    drhodpclim = (drho2./rg_dpr_clim + drho1./rg_dpr_clim_below)/2;
    pvclim = (1/rho_ref)*fcor.*drhodpclim; %pv climatological mean overall
    rg_pv = NaN(length(rg_xi),length(rg_yi),length(rg_yrgrid)); %calculate pv anomaly for R&G argo

    for it=1:length(rg_yrgrid)
        dpr2 = rg_dpr_anom(:,:,it); %anomaly of layer thickness between this layer and one above it
        dpr1 = rg_dpr_anom_below(:,:,it); %anomaly of layer thickness between this layer and one below it
        drhodp_tot = (drho2./(rg_dpr_clim(1:4:end,21:4:621) + dpr2) + drho1./(rg_dpr_clim_below(1:4:end,21:4:621) + dpr1))/2;
        rg_pv(:,:,it) = (1/rho_ref)*fcor.*drhodp_tot - pvclim(1:4:end,21:4:621);
    end 
    rg_pv = permute(rg_pv,[2,1,3]);

    rg_tseries = squeeze(rg_pv(idxx,idxy,is_rg,:)); 
    %now re-set the mean to zero now that I've restricted the
    %time series
    rg_seasonal_cycle = zeros(size(rg_tseries));
    
    %remove seasonal from R&G product
    %calculate seasonal cycle for this timeseries
    [zanom,zmean,coef,zm] = seasonal(rg_yrgrid,rg_tseries);
    rg_seasonal_cycle = zmean;
    rg_timeseries = rg_tseries - rg_seasonal_cycle;
    
    % CALCULATE ARGO PROFILE PV
    %use argo clim?
    
    is_dat = find(dat.sig1grid == sig1grid(is(ll)));
    
    %calculate profile data layer thickness from sigma1 grid, pressure; then calculate pv
    %use layer thickness climatology (pressure diff between sigma layer and the one above it)
    rg_dpr_clim = NaN(size(rg.pres_mean_sig1));
    for jj = 2:length(select_sig1grid)
        rg_dpr_clim(:,:,jj) =  rg.pres_mean_sig1(:,:,jj) - rg.pres_mean_sig1(:,:,jj-1);
    end
    rg_dpr_clim(:,:,1) = rg_dpr_clim(:,:,2); % just fill in the top sigma1 layer with layer below
    
    rg_dpr_anom = NaN(size(rg.pres_anom_sig1));
    for jj = 2:length(select_sig1grid) %for each sigma1 layer (36)
        for tt = 1:size(rg.pres_anom_sig1,4) %for each month
            rg_dpr_anom(:,:,jj,tt) = rg.pres_anom_sig1(:,:,jj,tt) - rg.pres_anom_sig1(:,:,jj-1,tt);
        end
    end
    rg_dpr_anom(:,:,1,:) = rg_dpr_anom(:,:,2,:); %just fill in the top sigma1 layer with layer below
    
    rg_dpr_clim_below = NaN(size(rg.pres_mean_sig1));
    for jj = 1:(length(select_sig1grid)-1)
        rg_dpr_clim_below(:,:,jj) = rg.pres_mean_sig1(:,:,jj+1) - rg.pres_mean_sig1(:,:,jj);
    end
    rg_dpr_clim_below(:,:,end) = rg_dpr_clim_below(:,:,2); % just fill in the bottom sigma1 layer with layer above
    
    rg_dpr_anom_below = NaN(size(rg.pres_anom_sig1));
    for jj = 1:(length(select_sig1grid)-1) %for each sigma1 layer (36)
        for tt = 1:size(rg.pres_anom_sig1,4) %for each month
            rg_dpr_anom_below(:,:,jj,tt) = rg.pres_anom_sig1(:,:,jj+1,tt) - rg.pres_anom_sig1(:,:,jj,tt);
        end
    end
    rg_dpr_anom_below(:,:,end,:) = rg_dpr_anom_below(:,:,end-1,:); %just fill in the bottom sigma1 layer with layer above
    
    fcor = calc_fcor(rg_xi,rg_yi);
    %drho2 and drho1 already calculated above for new argo product
    drhodpclim = (drho2./rg_dpr_clim + drho1./rg_dpr_clim_below)/2;
    pvclim = (1/rho_ref)*fcor.*drhodpclim; %pv climatological mean overall
    rg_pv = NaN(length(rg_xi),length(rg_yi),length(rg_yrgrid)); %calculate pv anomaly for R&G argo

    for it=1:length(rg_yrgrid)
        dpr2 = rg_dpr_anom(:,:,it); %anomaly of layer thickness between this layer and one above it
        dpr1 = rg_dpr_anom_below(:,:,it); %anomaly of layer thickness between this layer and one below it
        drhodp_tot = (drho2./(rg_dpr_clim(1:4:end,21:4:621) + dpr2) + drho1./(rg_dpr_clim_below(1:4:end,21:4:621) + dpr1))/2;
        rg_pv(:,:,it) = (1/rho_ref)*fcor.*drhodp_tot - pvclim(1:4:end,21:4:621);
    end 
    rg_pv = permute(rg_pv,[2,1,3]);

    rg_tseries = squeeze(rg_pv(idxx,idxy,is_rg,:)); 
    %now re-set the mean to zero now that I've restricted the
    %time series
    rg_seasonal_cycle = zeros(size(rg_tseries));
    
    %remove seasonal from R&G product
    %calculate seasonal cycle for this timeseries
    [zanom,zmean,coef,zm] = seasonal(rg_yrgrid,rg_tseries);
    rg_seasonal_cycle = zmean;
    rg_timeseries = rg_tseries - rg_seasonal_cycle;
    
    
    
    
    
    
    
    
    
    
    
    subplot(2,2,ll,'Position',[positions(ll,:) 0.43 0.33]);
    plot(dat.decyr(ii),dat.sa_sig(ii,is(ll)),'b.')
    hold on
    error_area(yrgrid,zz,3*ze,[1 1 1]*0.5) %stderr
    plot(yrgrid,zz,'k.-') %new loess
    plot(new_e_yrgrid,e_timeseries,':','Color',[0.62 0.03 0.92],'LineWidth',2); %ECCO v4
    plot(rg_yrgrid,rg_timeseries,'Color',[1 0.6 0.1],'LineWidth',1.5); %Roemmich & Gilson
    set(gca,'fontsize',11)
    xlim([2003 2024])
    if ll == 3 || ll == 4
        xlabel('year')
    end
    if ll == 1 || ll == 3
        ylabel('psu')
    end
    ylim([-0.5 0.5])
    title([name{ll},', sig1 =',num2str(sig1grid(is(ll)))],'FontSize',15)
   
end
leg1 = legend('profile data','95% confidence','loess','ECCO v4','R&G','Location','southoutside');
set(leg1,'Position',[0.4375 0.4175 0.1 0.1]);
sgtitle('PV Anomaly, seasonal cycle removed','FontSize',20)
saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/PV_timeseries_comparison.png');
saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/PV_timeseries_comparison.fig');

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

