%try making a scatter plot of advection time scale versus correlation 
%argo in one color, ecco in another

addpath /atalanta/home/swijffels/toolbox/susan
addpath /atalanta/home/swijffels/toolbox/eez
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/chersh/toolbox
cd /atalanta/home/chersh/SpiceAnomalies/persistence_maps
%% spice: corr vs time by sigma1 surface

filter_starting_variance = 0; %do you want to filter by starting variance
stdcutoff = 0.05;

filter_pressure = 0; %do you want to filter by pressure
pressure_range = [0 200]; %pressure range you want to use

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

for ss = 1:nsig
    figure(ss)
    clf
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    if filter_pressure == 1
        advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
        advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    end
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    plot([0 100],[0 0],'r','LineWidth',1)
    hold on
    s1 = scatter(advect_lag,advect_corr,5,[1 1 1]*0.4,'.');
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(lags,aa,bb,'k');
    alpha(ee,0.3)
    plot(lags,aa,'k','LineWidth',3)
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    if filter_pressure == 1
        advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
        advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    end
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    c2 = [95 2 235]./252;
    s2 = scatter(advect_lag,advect_corr,5,[167 132 219]./252,'.');
    hold on
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(lags,aa,bb,c2);
    alpha(ee,0.3)
    plot(lags,aa,'Color',c2,'LineWidth',3)
    xlim([0 96])
    ylim([-0.8 1.2])    
    title('spice: correlation vs time')
    
    if filter_starting_variance == 1
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_time_sig1_',sigma1,'_filter_variance.png']);
    elseif filter_pressure == 1
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_time_sig1_',sigma1,'_filter_pressure.png']);
    else
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_time_sig1_',sigma1,'.png']);
    end
end

%% spice: corr vs time by pressure range

filter_starting_variance = 0; %do you want to filter by starting variance
stdcutoff = 0.05;

pressure_range = [200 500]; %pressure range you want to use

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

space = 180000; %this is dumb but what is absolute max number of points each sig1 level could contribute
lagspace = 150000; %absolute max number of points each sig1 level could contribute to one lag
advect_corr_tot_e = NaN(space*nsig,1);
advect_lag_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(lagspace*nsig,length(lags));
advect_corr_tot_a = NaN(space*nsig,1);
advect_lag_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(lagspace*nsig,length(lags));
for ss = 1:nsig
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    advect_lag_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_lag;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(lagspace,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    advect_corr_matrix_tot_e(1+(ss-1)*lagspace:(ss*lagspace),:) = advect_corr_matrix;
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    advect_lag_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_lag;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    advect_corr_matrix_tot_a(1+(ss-1)*lagspace:(ss*lagspace),:) = advect_corr_matrix;
   
end

figure
%percentile of data you want to display
percent = [30,70];

plot([0 100],[0 0],'r','LineWidth',1)
hold on
s1 = scatter(advect_lag_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(lags,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(lags,aae,'k','LineWidth',3)
    
c2 = [95 2 235]./252;
s2 = scatter(advect_lag_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(lags,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(lags,aaa,'Color',c2,'LineWidth',3)
xlim([0 96])
ylim([-0.8 1.2])
ylabel('correlation')
xlabel('months')
title(['spice: pressure = ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])   

[h1,p1,ks2stat1] = kstest2(aaa,aae);

saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_time_filter_pressure',num2str(pressure_range(1)),'_',num2str(pressure_range(2)),'.png']);

%% spice: corr vs distance by sigma1 surface

filter_starting_variance = 0; %do you want to filter
stdcutoff = 0.05;

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

for ss = 1:nsig
    figure(ss)
    clf
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    %list of dist values in 100 km increments
    dists = 0:200:10^5; %should cut off really high distances??
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this lag value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
    end
    plot([0 15000],[0 0],'r','LineWidth',1)
    hold on
    s1 = scatter(distkm_round,advect_corr,5,[1 1 1]*0.4,'.');
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(dists,aa,bb,'k');
    alpha(ee,0.3)
    plot(dists,aa,'k','LineWidth',3)
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this lag value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
    end
    c2 = [95 2 235]./252;
    s2 = scatter(distkm_round,advect_corr,5,[167 132 219]./252,'.');
    hold on
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(dists,aa,bb,c2);
    alpha(ee,0.3)
    plot(dists,aa,'Color',c2,'LineWidth',3)
    xlim([0 15000])
    ylim([-0.8 1.2])
    title('spice: corr vs distance')

    [h1,p1,ks2stat1] = kstest2(aaa,aae);

    if filter_starting_variance == 1
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_dist_sig1_',sigma1,'_filter_variance.png']);
    else
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_dist_sig1_',sigma1,'.png']);
    end
end

%% spice: corr vs distance by pressure range

filter_starting_variance = 0; %do you want to filter
stdcutoff = 0.05;

pressure_range = [200 500]; %pressure range you want to use

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

%list of dist values in 200 km increments
dists = 0:200:10^5;

space = 180000; %this is dumb but what is absolute max number of points each sig1 level could contribute
distspace = 5000; %absolute max number of points each sig1 level could contribute to one distance bin
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

figure
%percentile of data you want to display
percent = [30,70];

plot([0 15000],[0 0],'r','LineWidth',1)
hold on
s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aae,'k','LineWidth',3)

c2 = [95 2 235]./252;
s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aaa,'Color',c2,'LineWidth',3)
xlim([0 15000])
ylim([-0.8 1.2])
title(['spice: corr vs distance, p = ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])   

[h1,p1,ks2stat1] = kstest2(aaa(5:end),aae(5:end));

saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/spice_corr_vs_dist_filter_pressure',num2str(pressure_range(1)),'_',num2str(pressure_range(2)),'.png']);

%% PV: corr vs time by sigma1 surface

filter_starting_variance = 0; %do you want to filter
stdcutoff = 4*10^(-11);

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('PV_argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('PV_ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

for ss = 1:nsig
    figure(ss)
    clf
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    plot([0 100],[0 0],'r','LineWidth',1)
    hold on
    s1 = scatter(advect_lag,advect_corr,5,[1 1 1]*0.4,'.');
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(lags,aa,bb,'k');
    alpha(ee,0.3)
    plot(lags,aa,'k','LineWidth',3)
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    c2 = [95 2 235]./252;
    s2 = scatter(advect_lag,advect_corr,5,[167 132 219]./252,'.');
    hold on
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(lags,aa,bb,c2);
    alpha(ee,0.3)
    plot(lags,aa,'Color',c2,'LineWidth',3)
    xlim([0 96])
    ylim([-0.8 1.2])
    xlabel('months')
    ylabel('correlation')
    title('PV: corr vs time')

    if filter_starting_variance == 1
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_time_sig1_',sigma1,'_filter_variance.png']);
    else
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_time_sig1_',sigma1,'.png']);
    end

end

%% PV: corr vs time by pressure range

filter_starting_variance = 0; %do you want to filter
stdcutoff = 3.5*10^(-11);

pressure_range = [0 1000]; %pressure range you want to use

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('PV_argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('PV_ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

space = 180000; %this is dumb but what is absolute max number of points each sig1 level could contribute
lagspace = 150000; %absolute max number of points each sig1 level could contribute to one lag
advect_corr_tot_e = NaN(space*nsig,1);
advect_lag_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(lagspace*nsig,length(lags));
advect_corr_tot_a = NaN(space*nsig,1);
advect_lag_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(lagspace*nsig,length(lags));
for ss = 1:nsig
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    advect_lag_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_lag;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(lagspace,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    advect_corr_matrix_tot_e(1+(ss-1)*lagspace:(ss*lagspace),:) = advect_corr_matrix;
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        advect_lag = advect_lag(sqrt(starting_variance) >= stdcutoff);
    end
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_lag = advect_lag(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    advect_lag_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_lag;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(lags));
    for ll = 1:length(lags)
    %all corr values that correspond to this lag value
        values = advect_corr(advect_lag == lags(ll));
        advect_corr_matrix(1:length(values),ll) = values;
    end
    advect_corr_matrix_tot_a(1+(ss-1)*lagspace:(ss*lagspace),:) = advect_corr_matrix;
   
end

figure
%percentile of data you want to display
percent = [30,70];

plot([0 100],[0 0],'r','LineWidth',1)
hold on
s1 = scatter(advect_lag_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aa = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(lags,aa,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(lags,aa,'k','LineWidth',3)
    
c2 = [95 2 235]./252;
s2 = scatter(advect_lag_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(lags,aa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(lags,aa,'Color',c2,'LineWidth',3)
xlim([0 96])
ylim([-0.8 1.2])
ylabel('correlation')
xlabel('months')
title(['PV: corr vs time, p = ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])   

saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_time_filter_pressure',num2str(pressure_range(1)),'_',num2str(pressure_range(2)),'.png']);

%% PV: corr vs distance by sigma1 surface

filter_starting_variance = 0; %do you want to filter
stdcutoff = 3.5*10^(-11);

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('PV_argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('PV_ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

test = 1;
for ss = 1:nsig
    figure(ss)
    clf
    load(e_sig1names{ss}) %plot ECCO
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    %list of dist values in 100 km increments
    dists = 0:200:10^5; %should cut off really high distances??
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this lag value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
    end
    plot([0 15000],[0 0],'r','LineWidth',1)
    hold on
    s1 = scatter(distkm_round,advect_corr,5,[1 1 1]*0.4,'.');
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(dists,aa,bb,'k');
    alpha(ee,0.3)
    plot(dists,aa,'k','LineWidth',3)
    
    load(a_sig1names{ss}) %plot Argo
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(150000,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this lag value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
    end
    c2 = [95 2 235]./252;
    s2 = scatter(distkm_round,advect_corr,5,[167 132 219]./252,'.');
    hold on
    aa = mean(advect_corr_matrix,1,'omitnan');
    bb = std(advect_corr_matrix,'omitnan');
    ee = error_area(dists,aa,bb,c2);
    alpha(ee,0.3)
    plot(dists,aa,'Color',c2,'LineWidth',3)
    xlim([0 15000])
    ylim([-0.8 1.2])
    xlabel('km')
    ylabel('correlation')
    
    title('PV: corr vs distance')

    if filter_starting_variance == 1
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_dist_sig1_',sigma1,'_filter_variance.png']);
    else
        saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_dist_sig1_',sigma1,'.png']);
    end
end

%% PV: corr vs distance by pressure range
filter_starting_variance = 0; %do you want to filter
stdcutoff = 3.5*10^(-11);

pressure_range = [0 1000]; %pressure range you want to use

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('PV_argo_down_streamline_persistence_sig*');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('PV_ecco_down_streamline_persistence_sig*');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

%list of dist values in 200 km increments
dists = 0:200:10^5;

space = 180000; %this is dumb but what is absolute max number of points each sig1 level could contribute
distspace = 5000; %absolute max number of points each sig1 level could contribute to one distance bin
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

figure
%percentile of data you want to display
percent = [30,70];

plot([0 15000],[0 0],'r','LineWidth',1)
hold on
s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aa = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aa,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aa,'k','LineWidth',3)

c2 = [95 2 235]./252;
s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aa,'Color',c2,'LineWidth',3)
xlim([0 15000])
ylim([-0.8 1.2])
title(['PV: corr vs distance, p = ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])   

saveas(gcf,['/atalanta/home/chersh/SpiceAnomalies/decorrelation_scatterplots/PV_corr_vs_dist_filter_pressure',num2str(pressure_range(1)),'_',num2str(pressure_range(2)),'.png']);

%% 4-panel plot: spice/PV on 0-200dbar and 200-500 dbar
cd /atalanta/home/chersh/SpiceAnomalies/persistence_maps
%do correlation vs distance
filter_starting_variance = 0; %do you want to filter
stdcutoff = 0.05;

figure
%plot spice first
addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('argo_down_streamline_persistence_sig*_tunnels.mat');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('ecco_down_streamline_persistence_sig*_tunnels.mat');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

maxlag = 96;
lags = -maxlag:maxlag;

%list of dist values in 200 km increments
dists = 0:200:10^5;

space = 180000; %this is dumb but what is absolute max number of points each sig1 level could contribute
distspace = 5000; %absolute max number of points each sig1 level could contribute to one distance bin

%percentile of data you want to display
percent = [30,70];
%-------------------------------------------------------------------------------
pressure_range = [0 200]; %pressure range you want to use
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

subaxis(2,2,1,'SV',0.03,'SH',0.03)

plot([0 10000],[0 0],'r','LineWidth',1)
hold on
%s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aae,'k','LineWidth',3)

c2 = [95 2 235]./252;
%s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aaa,'Color',c2,'LineWidth',3)
xlim([0 10000])
ylim([-0.2 1.2])
ylabel('correlation')
set(gca,'XTick',[])
vax = axis;
text(vax(1)+2000,vax(4)-0.2,['Spice: ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])

%compare argo and ecco with kstest
h1 = NaN(size(dists));
p1 = NaN(size(dists));
ks2stat1 = NaN(size(dists));
for dd = 1:length(dists)
    if sum(~isnan(advect_corr_matrix_tot_e(:,dd))) > 5 && sum(~isnan(advect_corr_matrix_tot_a(:,dd))) > 5
        [h1(dd),p1(dd),ks2stat1(dd)] = kstest2(advect_corr_matrix_tot_a(:,dd),advect_corr_matrix_tot_e(:,dd));
    end
end

%-------------------------------------------------------------------------------
pressure_range = [200 500]; %pressure range you want to use
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr > pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr > pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

subaxis(2,2,3,'SV',0.03,'SH',0.03)

plot([0 10000],[0 0],'r','LineWidth',1)
hold on
%s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aae,'k','LineWidth',3)

c2 = [95 2 235]./252;
%s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aaa,'Color',c2,'LineWidth',3)
xlim([0 10000])
ylim([-0.2 1.2])
ylabel('correlation')
xlabel('Distance from outcrop (km)')
vax = axis;
text(vax(1)+2000,vax(4)-0.2,['Spice: ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])

%compare argo and ecco with kstest
h2 = NaN(size(dists));
p2 = NaN(size(dists));
ks2stat2 = NaN(size(dists));
for dd = 1:length(dists)
    if sum(~isnan(advect_corr_matrix_tot_e(:,dd))) > 5 && sum(~isnan(advect_corr_matrix_tot_a(:,dd))) > 5
        [h2(dd),p2(dd),ks2stat2(dd)] = kstest2(advect_corr_matrix_tot_a(:,dd),advect_corr_matrix_tot_e(:,dd));
    end
end

%now plot PV
%------------------------------------------------------------------------------
filter_starting_variance = 0; %do you want to filter
stdcutoff = 3.5*10^(-11);

addpath('/atalanta/home/chersh/SpiceAnomalies/persistence_maps');
a_sig1files = dir('PV_argo_down_streamline_persistence_sig*_tunnels.mat');
nsig = length(a_sig1files); %number of years of data
a_sig1names = cell(size(a_sig1files));
for ii = 1:length(a_sig1files) %probably a better way to do this but w/e
    a_sig1names{ii} = a_sig1files(ii).name;
end

e_sig1files = dir('PV_ecco_down_streamline_persistence_sig*_tunnels.mat');
%nsig = length(e_sig1files); %number of years of data
e_sig1names = cell(size(e_sig1files));
for ii = 1:length(e_sig1files) %probably a better way to do this but w/e
    e_sig1names{ii} = e_sig1files(ii).name;
end

pressure_range = [0 200]; %pressure range you want to use
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

subaxis(2,2,2,'SH',0.03,'SV',0.03)

plot([0 10000],[0 0],'r','LineWidth',1)
hold on
%s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aae,'k','LineWidth',3)

c2 = [95 2 235]./252;
%s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aaa,'Color',c2,'LineWidth',3)
xlim([0 10000])
ylim([-0.2 1.2])
set(gca,'YTick',[])
set(gca,'XTick',[])
vax = axis;
text(vax(1)+2000,vax(4)-0.2,['PV: ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])

h3 = NaN(size(dists));
p3 = NaN(size(dists));
ks2stat3 = NaN(size(dists));
for dd = 1:length(dists)
    if sum(~isnan(advect_corr_matrix_tot_e(:,dd))) > 5 && sum(~isnan(advect_corr_matrix_tot_a(:,dd))) > 5
        [h3(dd),p3(dd),ks2stat3(dd)] = kstest2(advect_corr_matrix_tot_a(:,dd),advect_corr_matrix_tot_e(:,dd));
    end
end
%-------------------------------------------------------------------------------------
pressure_range = [200 500]; %pressure range you want to use
advect_corr_tot_e = NaN(space*nsig,1);
distkm_round_tot_e = NaN(space*nsig,1);
advect_corr_matrix_tot_e = NaN(distspace*nsig,length(dists));
advect_corr_tot_a = NaN(space*nsig,1);
distkm_round_tot_a = NaN(space*nsig,1);
advect_corr_matrix_tot_a = NaN(distspace*nsig,length(dists));
for ss = 1:nsig
    %ECCO
    load(e_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr > pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr > pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_e(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;

    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    test = 1;
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_e(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;
   
    %ARGO
    load(a_sig1names{ss})
    if filter_starting_variance == 1
        advect_corr = advect_corr(sqrt(starting_variance) >= stdcutoff);
        distkm = distkm(sqrt(starting_variance) >= stdcutoff);
    end
    %round distkm to nearest 100 km
    distkm_round = round(distkm,-2);
    
    %filter by pressure range
    advect_corr = advect_corr(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    distkm_round = distkm_round(mean_pr >= pressure_range(1) & mean_pr <= pressure_range(2));
    advect_corr_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = advect_corr;
    distkm_round_tot_a(1+(ss-1)*space:(ss-1)*space+length(advect_corr)) = distkm_round;
    
    %organize correlation values into a matrix of data
    advect_corr_matrix = NaN(distspace,length(dists));
    for dd = 1:length(dists)
    %all corr values that correspond to this distance value
        values = advect_corr(distkm_round == dists(dd));
        advect_corr_matrix(1:length(values),dd) = values;
        if length(values) > test
            test = length(values);
        end
    end
    advect_corr_matrix_tot_a(1+(ss-1)*distspace:(ss*distspace),:) = advect_corr_matrix;

end

subaxis(2,2,4,'SH',0.03,'SV',0.03)

plot([0 10000],[0 0],'r','LineWidth',1)
hold on
%s1 = scatter(distkm_round_tot_e,advect_corr_tot_e,5,[1 1 1]*0.4,'.');
aae = median(advect_corr_matrix_tot_e,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_e,percent,1);
ee = percentile_area(dists,aae,bb(1,:),bb(2,:),'k');
alpha(ee,0.3)
plot(dists,aae,'k','LineWidth',3)

c2 = [95 2 235]./252;
%s2 = scatter(distkm_round_tot_a,advect_corr_tot_a,5,[167 132 219]./252,'.');
aaa = median(advect_corr_matrix_tot_a,1,'omitnan');
bb = prctile(advect_corr_matrix_tot_a,percent,1);
ee = percentile_area(dists,aaa,bb(1,:),bb(2,:),c2);
alpha(ee,0.3)
plot(dists,aaa,'Color',c2,'LineWidth',3)
xlim([0 10000])
ylim([-0.2 1.2])
set(gca,'YTick',[])
xlabel('Distance from outcrop (km)')
vax = axis;
text(vax(1)+2000,vax(4)-0.2,['PV: ',num2str(pressure_range(1)),' to ',num2str(pressure_range(2)),' dbar'])
legend('','ECCO 30th-70th prctile','ECCO median','Argo 30th-70th prctile','Argo median')

h4 = NaN(size(dists));
p4 = NaN(size(dists));
ks2stat4 = NaN(size(dists));
for dd = 1:length(dists)
    if sum(~isnan(advect_corr_matrix_tot_e(:,dd))) > 5 && sum(~isnan(advect_corr_matrix_tot_a(:,dd))) > 5
        [h4(dd),p4(dd),ks2stat4(dd)] = kstest2(advect_corr_matrix_tot_a(:,dd),advect_corr_matrix_tot_e(:,dd));
    end
end
