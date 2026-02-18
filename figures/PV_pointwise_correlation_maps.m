%copied over from atalanta/home/chersh/SpiceAnomalies/argo_ecco_pointwise_correlation.m
%% add necessary paths
addpath /atalanta/home/swijffels/work/argo/
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/toolbox
addpath /atalanta/home/chersh/toolbox/cmocean
cd /atalanta/home/swijffels/work/argo/gridNSF

%% find data files
%load in Argo anomalies
if ~exist('a_despike')
    a_despike = matfile('/atalanta/home/chersh/SpiceAnomalies/despike_argo_data.mat');
end

%load in ECCO output
if ~exist('e')
    e = matfile('/batou/ECCOv4r4/exps/iter129_bulkformula/run/regularpoles/mat_current/sig1_gridded/gridonSig1_ecco4r4_iter129_bulkformula.mat');   
end

sig1grid = a_despike.sig1grid; %sigma1 level grid
a_yrgrid = a_despike.yrgrid; %months
a_xi = a_despike.xi; %longitude grid
a_yi = a_despike.yi; %latitude grid

e_sig1grid = e.sig1grid;
e_yrgrid = e.yrgrid;
e_xi = e.xi;
e_yi = e.yi;

load('/atalanta/home/chersh/SpiceAnomalies/argo_ecco_pointwise_correlation/argo_ecco_correlation_rem_linear_seasonal_Jun_2023.mat')

%% plot (maps on sig1 levels)
addpath /atalanta/home/chersh/toolbox/cmocean
addpath /atalanta/home/chersh/toolbox/stipple

figure('Position',[10 10 1200 600])
ax = 1;
cmap = cmocean('balance',12);
cmap = cmap(2:end-1,:);
for isig = [28,38,44,54]
    disp(ax)
    %layer thickness climatology (pressure diff between sigma layer and the one above it)
    dpclim2 = squeeze(clim.dpr.m(:,:,isig)) + squeeze(clim.dpr.c(:,:,isig,1));
    %layer thickness climatology (pressure diff between sigma layer and the one below it)
    dpclim1 = squeeze(clim.dpr.m(:,:,isig+1)) + squeeze(clim.dpr.c(:,:,isig+1,1));
    rho_ref = 1000;
    fcor = calc_fcor(size(dpclim2,1),size(dpclim2,2)); %coriolis parameter
    
    drho2 = sig1grid(isig) - sig1grid(isig-1); %change in density between this layer and the one above it
    drho1 = sig1grid(isig+1) - sig1grid(isig); %change in density between this layer and the one below it

    drhodpclim = (drho2./dpclim2 + drho1./dpclim1)/2;
    mean_field = (1/rho_ref)*fcor.*drhodpclim; %pv climatological mean overall
    
    lims = [30,70];
    data = squeeze(argo_ecco_corrcoef_pv(:,:,isig));
    interval = 0.1;
    mask = squeeze(pvalue_pv(:,:,isig)) <= 0.05;

    subaxis(2,2,ax,'SpacingVert',0.01,'SpacingHoriz',0.01)
    levels = -1:0.2:1;
    contourf(a_xi,a_yi,data,levels,'edgecolor','none')
    hold on
    
    [XI,YI] = meshgrid(a_xi,a_yi);
    stipple(XI,YI,mask,'density',200,'markersize',2);
    
    gebco('k');
    
    yticks([-50 -40 -30 -20 -10 0 10 20 30 40 50])
    yticklabels({'50S','40S','30S','20S','10S','0','10N','20N','30N','40N','50N'})
    xticklabels({'0','50E','100E','150E','160W','110W','60W','10W'})
    
    if ax == 2
       c1 = colorbar;
       c1.Position = [0.91 0.1 0.02 0.8];
       c1.Orientation = 'Vertical';
    end
    if ax == 1 || ax == 2
        set(gca,'XTick',[])
    end
    if ax == 2 || ax == 4
        set(gca,'YTick',[])
    end
    colormap(cmap)
    caxis([levels(1) levels(end)]);
    for xx = 1:10:length(a_xi)
        %plot([a_xi(xx),a_xi(xx)],[a_yi(1),a_yi(end)],'k:');
    end
    for yy = 6:10:length(a_yi)
        %plot([a_xi(1),a_xi(end)],[a_yi(yy),a_yi(yy)],'k:');
    end
    ylim([-55 55])
    vax = axis;
    text(vax(1)+57,vax(4)-22,['\sigma_1 = ',num2str(sig1grid(isig))],'FontSize',8)
    ax = ax + 1;
end

saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/pv_pointwise_correlation.png');
saveas(gcf,'/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/pv_pointwise_correlation.fig');

%% make some useful functions so this code isn't so long...

function [sig1grid,a_yrgrid,a_xi,a_yi,e_yrgrid,e_xi,e_yi,e_tidx,a_tidx] = startup(a,e)
    addpath /atalanta/home/swijffels/work/argo/

    addpath /atalanta/home/swijffels/toolbox/seawater
    addpath /atalanta/home/swijffels/toolbox/csirolib
    addpath /atalanta/home/swijffels/work/argo/matlab
    addpath /atalanta/home/swijffels/toolbox/subaxis
    addpath /atalanta/home/chersh/toolbox
    cd /atalanta/home/swijffels/work/argo/gridNSF

    %load in Argo anomalies
    %will have to load in Argo and ECCO data layer by layer and save, then
    %clear variables, because I can't load in everything at the same time

    sig1grid = a.sig1grid; %sigma1 level grid
    a_yrgrid = a.yrgrid; %months
    a_xi = a.xi; %longitude grid
    a_yi = a.yi; %latitude grid

    %load in ECCO output
    e_yrgrid = e.yrgrid;
    e_xi = e.xi;
    e_yi = e.yi;

    %find overlap in time between Argo and ECCO
    %Argo start is start bound, ECCO end is end bound
    e_tidx = e_yrgrid > a_yrgrid(1);
    a_tidx = a_yrgrid < e_yrgrid(end);

end 

function a_var = load_a_anom(a,var,isig) %this just loads argo anom data for 'sa', 'pr', or 'dpr'

    if strcmp(var,'sa') == 1
        a_var = squeeze(a.sa_fwa(:,:,isig,:));
        a_var_wa = squeeze(a.sa_wa(:,:,isig,:));  
    elseif strcmp(var,'pr') == 1
        a_var = squeeze(a.pr_fwa(:,:,isig,:));
        a_var_wa = squeeze(a.pr_wa(:,:,isig,:));  
    elseif strcmp(var,'dpr') == 1
        a_var = squeeze(a.dpr_fwa(:,:,isig,:));
        a_var_wa = squeeze(a.dpr_wa(:,:,isig,:));  
    end
    a_ir = isnan(a_var) & ~isnan(a_var_wa);
    a_var(a_ir) = a_var_wa(a_ir);
    clear *_wa *_fwa 
end

function fcor = calc_fcor(xgrid,ygrid)
    fcor = NaN(length(xgrid),length(ygrid)); %coriolis parameter
    for aa = 1:length(xgrid)
        for bb = 1:length(ygrid)
            lati = ygrid(bb);
            fcor(aa,bb) = 2*7.2921*10^(-5)*sind(lati); %fill out coriolis parameter array
        end
    end
end

function fill_a_var = fill_a(a_var,a_xi,a_yi,a_yrgrid) 
    location_mask = logical.empty(length(a_xi),length(a_yi),0);

    for xx = 1:length(a_xi)
        for yy = 1:length(a_yi)
            if sum(isnan(a_var(xx,yy,:)))/length(a_yrgrid) < 0.50 % mask out if 50% or more of data at this location is blank
                location_mask(xx,yy,1) = false;
            else
                location_mask(xx,yy,1) = true;
            end
        end
    end
    
    %use inpaint_nans function from matlab file exchange to do 2D spatial
    %interpolation of Argo data for each time step
    fill_a_var = NaN(length(a_xi),length(a_yi),length(a_yrgrid));
    for tt = 1:length(a_yrgrid)
        fill_a_var(:,:,tt) = inpaint_nans(a_var(:,:,tt));
    end
    %now apply the mask
    for tt = 1:length(a_yrgrid)
        junk = fill_a_var(:,:,tt); %this seems dumb but indexing not working otherwise
        junk(location_mask) = NaN;
        fill_a_var(:,:,tt) = junk;
    end
    
    fill_a_var = permute(fill_a_var,[2 1 3]); %just transposing this to match ecco 

end

function [a_timeseries,new_a_yrgrid] = get_a_timeseries(tseries,a_yrgrid,e_yrgrid,rem_seas,rem_lin) %pulls out argo timeseries, restricts to overlap w ECCO (and de-means), removes linear and/or seasonal variation
    a_tseries = squeeze(tseries);
    new_a_yrgrid = a_yrgrid(~isnan(a_tseries)); %cut off yrgrid corresponding to any NaNs at the beginning
    a_tseries = a_tseries(~isnan(a_tseries)); %cut off data NaNs at the beginning
    a_tidx = new_a_yrgrid <= e_yrgrid(end) & new_a_yrgrid >= 2004; %cut off before 2004 if not done already
    a_tseries = a_tseries(a_tidx); %restrict data to overlap with ecco
    new_a_yrgrid = new_a_yrgrid(a_tidx); % restrict yrgrid to overlap with ecco
    %now re-set the mean to zero now that I've restricted the
    %time series
    a_tseries = a_tseries - mean(a_tseries,'omitnan');
    a_seasonal_cycle = zeros(size(a_tseries));
    a_linear_trend = zeros(size(a_tseries));
                
    if rem_seas == 1
    %remove seasonal from Argo  
        %calculate seasonal cycle for this timeseries
        [~,zmean,~,~] = seasonal(new_a_yrgrid,a_tseries);
        a_seasonal_cycle = zmean;
    end
            
    if rem_lin == 1
        %calculate linear trend for Argo
        %need to interpolate or give query points in the case of NaNs
        p = polyfit(new_a_yrgrid,a_tseries,1);
        a_linear_trend = polyval(p,new_a_yrgrid);
        
    end
    
    a_timeseries = a_tseries - a_seasonal_cycle - a_linear_trend';

end

function e_timeseries = get_e_timeseries(tseries,e_yrgrid,new_a_yrgrid,rem_seas,rem_lin) %pulls out ECCO timeseries, restricts to overlap w Argo (and de-means), removes linear and/or seasonal variation
    
    e_tseries = squeeze(tseries);
    e_tseries = interp1(e_yrgrid,e_tseries,new_a_yrgrid);  %first interpolate ECCO timeseries to same time of month as Argo data
    %e_tidx = e_yrgrid >= new_a_yrgrid(1) & e_yrgrid <= new_a_yrgrid(end); %cut off ecco before Argo data starts (or after Argo data ends in a couple of cases)
    
    %new_e_yrgrid = e_yrgrid(e_tidx); %restrict yrgrid to overlap with argo
    %e_tseries = e_tseries(e_tidx); %restrict data to overlap with argo
    
    %now re-set the mean to zero now that I've restricted the
    %time series
    e_tseries = e_tseries - mean(e_tseries,'omitnan');
    e_seasonal_cycle = zeros(size(e_tseries));
    e_linear_trend = zeros(size(e_tseries));
    
    if rem_seas == 1
    %remove seasonal from ECCO  
        %calculate seasonal cycle for this timeseries
        [~,zmean,~,~] = seasonal(new_a_yrgrid,e_tseries);
        e_seasonal_cycle = zmean;
    end
            
    if rem_lin == 1
        %calculate linear trend for ECCO
        %need to interpolate or give query points in the case of NaNs
        yrgrid_nonan = new_a_yrgrid(~isnan(e_tseries));
        e_tseries_nonan = e_tseries(~isnan(e_tseries));
        fill_tseries = interp1(yrgrid_nonan,e_tseries_nonan,new_a_yrgrid); %check this
        fill_tseries(end) = fill_tseries(end-1); %somehow the last one was always NaN
        p = polyfit(new_a_yrgrid,fill_tseries,1);
        e_linear_trend = polyval(p,new_a_yrgrid);
        %calculate linear trend for ECCO
        %need to interpolate or give query points in the case of NaNs
    end
    
    e_timeseries = e_tseries' - e_seasonal_cycle - e_linear_trend';

end

function pvalue = montecarlo_significance(argo_timeseries,ecco_timeseries)
    % calculate degree of freedom? class notes eqn 2.216
    %follow procedure from pset 4 of 12.805, question 2 for a few locations to see how
    %pvalue changes from the corrcoef value
    
    %what is the correlation and pvalue (assuming monthly independence) of this
    %pair of time series?
    [Corr,Pval] = corrcoef(argo_timeseries,ecco_timeseries);
    Corr = Corr(1,2);
    Pval = Pval(1,2);

    %Cholesky decomposition of spice and PV time series
    [C_argo,~] = xcov(argo_timeseries,'coeff');
    C_argo = C_argo(length(argo_timeseries):end); %maybe I don't want the negative lags?
    [C_ecco,~] = xcov(ecco_timeseries,'coeff');
    C_ecco = C_ecco(length(ecco_timeseries):end); %maybe I don't want the negative lags?
    Cxx_argo = toeplitz(C_argo);
    Cxx_ecco = toeplitz(C_ecco);
    small = 10^-15;
    test = 0;
    try chol(Cxx_argo + small*eye(size(Cxx_argo)));
        %disp('symmetric pos def')
    catch ME
        disp('not symmetric pos def')
        test = 1;
    end
    if test == 0
        R_argo = chol(Cxx_argo + small*eye(size(Cxx_argo)));
    end
    
    test = 0;
    try chol(Cxx_ecco + small*eye(size(Cxx_ecco)));
        %disp('symmetric pos def')
    catch ME
        disp('not symmetric pos def')
        test = 1;
    end
    if test == 0
        R_ecco = chol(Cxx_ecco + small*eye(size(Cxx_ecco)));
    end

    %run 1k spice/pv pairs with this autocovariance and calculate their correlations
    numtests = 1000;
    correlations = NaN(numtests,1);
    for tt = 1:length(correlations)
        rand_argo = R_argo'*randn(length(R_argo),1);
        rand_ecco = R_ecco'*randn(length(R_ecco),1);
        CC = corrcoef(rand_argo,rand_ecco);
        correlations(tt) = CC(1,2);
    end

    [counts,~] = histcounts(correlations);
    cdf = cumsum(counts/numtests);
    xaxis = linspace(-1,1,length(cdf));
    [~,corrlocation] = min(abs(xaxis-Corr));
    if Corr > 0
        newpval = 2*(1-cdf(corrlocation)); %find p value from the cdf based on actual correlation of spice and pv in this location
    else
        newpval = 2*(cdf(corrlocation));
    end
    
    pvalue = newpval;
end


