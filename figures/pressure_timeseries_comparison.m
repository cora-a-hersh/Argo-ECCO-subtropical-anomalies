addpath /atalanta/home/swijffels/toolbox/susan
addpath /atalanta/home/swijffels/toolbox/eez
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
cd /atalanta/home/swijffels/work/argo/gridNSF

%% load in data
if ~exist('sa_fwa')
anom = matfile('gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
end
yrgrid = anom.yrgrid;
a_xi = anom.xi;
a_yi = anom.yi;
sig1grid = anom.sig1grid;

%% get input data
dat = load('Sigma1_anomalies_Argo_CORA_sla_superobs_1998_2020.mat','decyr','sa_sig','sas_sig','lat','lon','pa_sig','pas_sig');

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

%% time series pressure
addpath /atalanta/home/chersh/toolbox

xplot = [220,320,220,80];
yplot = [20,20,-20,-20];
is = [38,44,34,40]; %indices of sigma1 levels
name = {'North Pacific','North Atlantic','South Pacific','South Indian'};

figure('Position',[10 10 1600 900])

for ll = 1:length(xplot)
    ii = abs(dat.lon - xplot(ll))<2 & abs(dat.lat - yplot(ll))<1.5;
    
    listx = 1:length(a_xi); %this is stupid
    listy = 1:length(a_yi);
    idxx = listx(a_xi==xplot(ll));
    idxy = listy(a_yi==yplot(ll));
    zz = squeeze(anom.pr_fwa(idxx,idxy,is(ll),:));
    ze = squeeze(anom.pr_fwe(idxx,idxy,is(ll),:));

    za = squeeze(anom.pr_wa(idxx,idxy,is(ll),:));
    
    % get ecco timeseries
    [~,idxx] = min(abs(e_xi - xplot(ll)));
    [~,idxy] = min(abs(e_yi - yplot(ll)));
    e_tseries = squeeze(e.pri(idxx,idxy,is(ll)-13,:));
    e_tseries = e_tseries(e_tidx); %restrict to overlap with argo
    %now re-set the mean to zero now that I've restricted the
    %time series
    e_tseries = e_tseries - mean(e_tseries,'omitnan');
    e_seasonal_cycle = zeros(size(e_tseries));
    new_e_yrgrid = e_yrgrid(e_tidx); % restrict yrgrid to overlap with ecco
    
    %remove seasonal from ECCO  
    %calculate seasonal cycle for this timeseries
    [zanom,zmean,coef,zm] = seasonal(new_e_yrgrid,e_tseries);
    e_seasonal_cycle = zmean;
    e_timeseries = e_tseries - e_seasonal_cycle;
    
    % get R&G timeseries
    [~,idxx] = min(abs(rg_xi - xplot(ll)));
    [~,idxy] = min(abs(rg_yi - yplot(ll)));
    is_rg = find(select_sig1grid == sig1grid(is(ll)));
    rg_tseries = squeeze(rg.pres_anom_sig1(idxx,idxy,is_rg,:));
    %now re-set the mean to zero now that I've restricted the
    %time series
    rg_seasonal_cycle = zeros(size(rg_tseries));
    
    %remove seasonal from R&G product
    %calculate seasonal cycle for this timeseries
    [zanom,zmean,coef,zm] = seasonal(rg_yrgrid,rg_tseries);
    rg_seasonal_cycle = zmean;
    rg_timeseries = rg_tseries - rg_seasonal_cycle;
    
    subaxis(2,2,ll,'SpacingHoriz',0.035,'SpacingVert',0.06);
    plot(dat.decyr(ii),dat.pa_sig(ii,is(ll)),'b.','LineWidth',1)
    hold on
    plot(rg_yrgrid,rg_timeseries,'Color',[0.28 0.72 0.40],'LineWidth',1.5); %Roemmich & Gilson
    error_area(yrgrid,zz,3*ze,[1 1 1]*0.5)
    plot(yrgrid,zz,'k','LineWidth',2)
    plot(new_e_yrgrid,e_timeseries,'Color',[1 0.6 0.1],'LineWidth',2);
    set(gca,'fontsize',11)
    xlim([2004 2022])
    ylim([-65 65])
    if ll == 3 || ll == 4
        xlabel('year')
    end
    if ll == 1 || ll == 3
        ylabel('dbar')
        yticks([-60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60])
    end
    if ll == 1 || ll == 2
        set(gca,'XTicklabel',[]);
    end
    if ll == 2 || ll == 4
        set(gca,'YTicklabel',[]);
    end
    title(name{ll},'FontSize',15)
    vax = axis;
    text(vax(1)+1,vax(4)-10,['\sigma_1 = ',num2str(sig1grid(is(ll)))],'FontSize',12)
end
leg1 = legend('profile data','R&G','95% confidence','ArgoLoS','ECCO','Location','southoutside');
set(leg1,'Position',[0.4375 0.4175 0.1 0.1]);

sgtitle('Pressure anomaly, seasonal cycle removed','FontSize',20)
saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/pressure_timeseries_comparison.png']);
saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/Paper1_figures/pressure_timeseries_comparison.fig']);
