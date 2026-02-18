%put R&G Argo product on sigma1 grid for comparison to ECCO, other Argo
%products
%only do when actually on atalanta, because Ben's home folder isn't mounted
%on taris

addpath /home/bgreenwood/website/gridded
addpath /home/swijffels/toolbox/susan
addpath /home/swijffels/toolbox/eez
addpath /home/swijffels/toolbox/seawater
addpath /home/swijffels/toolbox/csirolib
addpath /home/swijffels/work/argo/matlab
addpath /home/swijffels/toolbox/subaxis
cd /home/swijffels/work/argo/gridNSF

anom = matfile('gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
sig1grid = anom.sig1grid;

%pick a few sigma1 level indices to save
sigma1_select = 21:56;
%sigma1_select = [21,28,34,35,38,40,41,44,48,54,56];
select_sig1grid = sig1grid(sigma1_select);

psal_anom = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','ARGO_SALINITY_ANOMALY');
psal_mean = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','ARGO_SALINITY_MEAN'); %mean from 2004-2018
temp_anom = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Temp.nc','ARGO_TEMPERATURE_ANOMALY');
temp_mean = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Temp.nc','ARGO_TEMPERATURE_MEAN'); %mean from 2004-2018

RG_xi = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','LONGITUDE');
RG_yi = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','LATITUDE');
RG_pri = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','PRESSURE');
RG_pri = double(RG_pri);
RG_time = ncread('/home/bgreenwood/website/gridded/RG_ArgoClim_Psal.nc','TIME'); %months since 2004-01-01 00:00:00

%for each vertical profile, calculate potential density referenced to 1000
%m using seawater toolbox
psal_anom_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select),length(RG_time));
psal_mean_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select));
temp_anom_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select),length(RG_time));
temp_mean_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select));
pres_anom_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select),length(RG_time));
pres_mean_sig1 = NaN(length(RG_xi),length(RG_yi),length(sigma1_select));
for xx = 1:length(RG_xi) %lon
    disp(xx)
    for yy = 1:length(RG_yi) %lat
        psal_mean_profile = squeeze(psal_mean(xx,yy,:));
        temp_mean_profile = squeeze(temp_mean(xx,yy,:));
        sigma1_profile = sw_pden(psal_mean_profile,temp_mean_profile,RG_pri,1000);
        if all(isnan(sigma1_profile))
            continue
        end
        sigma1_profile = sigma1_profile - 1000; % subtract 1000 from sigma1
        psal_mean_profile = psal_mean_profile(~isnan(sigma1_profile));
        temp_mean_profile = temp_mean_profile(~isnan(sigma1_profile));
        pres_mean_profile = RG_pri(~isnan(sigma1_profile));
        sigma1_profile = sigma1_profile(~isnan(sigma1_profile));
        psal_mean_sig1(xx,yy,:) = interp1(sigma1_profile,psal_mean_profile,select_sig1grid);
        temp_mean_sig1(xx,yy,:) = interp1(sigma1_profile,temp_mean_profile,select_sig1grid);
        pres_mean_sig1(xx,yy,:) = interp1(sigma1_profile,pres_mean_profile,select_sig1grid);
        for tt = 1:length(RG_time) %months
            psal_profile = squeeze(psal_anom(xx,yy,:,tt)) + squeeze(psal_mean(xx,yy,:)); %total salinity
            temp_profile = squeeze(temp_anom(xx,yy,:,tt)) + squeeze(temp_mean(xx,yy,:)); %total temp
            sigma1_profile = sw_pden(psal_profile,temp_profile,RG_pri,1000);
            if all(isnan(sigma1_profile))
                continue
            end
            sigma1_profile = sigma1_profile - 1000;
            psal_anom_profile = squeeze(psal_anom(xx,yy,:,tt));
            psal_anom_profile = psal_anom_profile(~isnan(sigma1_profile)); %get rid of NaNs
            temp_anom_profile = squeeze(temp_anom(xx,yy,:,tt));
            temp_anom_profile = temp_anom_profile(~isnan(sigma1_profile)); %get rid of NaNs
            pres_anom_profile = RG_pri(~isnan(sigma1_profile)); %get rid of NaNs in pressure profile
            sigma1_profile = sigma1_profile(~isnan(sigma1_profile)); %get rid of NaNs
            psal_anom_sig1(xx,yy,:,tt) = interp1(sigma1_profile,psal_anom_profile,select_sig1grid);
            temp_anom_sig1(xx,yy,:,tt) = interp1(sigma1_profile,temp_anom_profile,select_sig1grid);
            pres_anom_sig1(xx,yy,:,tt) = interp1(sigma1_profile,pres_anom_profile,select_sig1grid);
        end
    end
end

xi = RG_xi;
yi = RG_yi;
yrgrid = 2004 + RG_time*(1/12);
save('/home/chersh/SpiceAnomalies/RandG_on_sigma1.mat','psal_mean_sig1','temp_mean_sig1','pres_mean_sig1','psal_anom_sig1','temp_anom_sig1','pres_anom_sig1','select_sig1grid','xi','yi','yrgrid','-v7.3');

