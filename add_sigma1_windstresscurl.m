%% calculate sigma1, add to P gridded file

expname = 'nointerannual_wind';
%rootdir = '/vast'; %on taris or ilko; or /batou for iter129_bulkformula, nointerannual, and interannual_southpac
rootdir = '/vast/proj/ecco'; %vast path is different on atalanta
%rootdir = '/batou';
%addpath /atalanta/home/swijffels/toolbox/seawater %on taris or ilko
addpath /home/swijffels/toolbox/seawater %on atalanta

e = matfile([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/mat/P_gridded/gridonP_ecco4r4_',expname,'.mat']);
xi = e.xi;
yi = e.yi;
zi = e.zi;
yrgrid = e.yrgrid;
%%
%{
sig1i = NaN(size(e.pti));
for ii = 1:length(e.yrgrid) 
    disp(ii)
    for yy = 1:length(e.yi)
        si = squeeze(e.si(:,yy,:,ii));
        pti = squeeze(e.pti(:,yy,:,ii));
        ti = sw_temp(si,pti,zi',0);
        sig1i(:,yy,:,ii) = sw_pden(si,ti,zi',1000);
    end
end
%}

%hold salinity to background mean, let temperature vary --> find surface
%sigma1
sig1i_sconst = NaN(length(xi),length(yi),length(yrgrid));
si = squeeze(e.si(:,:,1,:));
pti = squeeze(e.pti(:,:,1,:));
si_mean = mean(si,3,'omitnan');
for ii = 1:length(yrgrid) 
    disp(ii)
    for yy = 1:length(yi)
        S = squeeze(si_mean(:,yy));
        PT = squeeze(pti(:,yy,ii));
        T = sw_temp(S,PT,zi(1),0);
        sig1i_sconst(:,yy,ii) = sw_pden(S,T,zi(1),1000);
    end
end

%hold temperature  to background mean, let salinity vary --> find surface
%sigma1
sig1i_tconst = NaN(length(xi),length(yi),length(yrgrid));
si = squeeze(e.si(:,:,1,:));
pti = squeeze(e.pti(:,:,1,:));
pti_mean = mean(pti,3,'omitnan');
for ii = 1:length(yrgrid) 
    disp(ii)
    for yy = 1:length(yi)
        S = squeeze(si(:,yy,ii));
        PT = squeeze(pti_mean(:,yy));
        T = sw_temp(S,PT,zi(1),0);
        sig1i_tconst(:,yy,ii) = sw_pden(S,T,zi(1),1000);
    end
end

save([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/mat/P_gridded/gridonP_ecco4r4_',expname,'.mat'],'sig1i_tconst','sig1i_sconst','-append');

%%
sconst = e.sig1i_sconst;
tconst = e.sig1i_tconst;
tot = squeeze(e.sig1i(:,:,1,:));

%%
sconst_sep = sconst(:,:,9:12:end);
sconst_mar = sconst(:,:,3:12:end);
tconst_sep = tconst(:,:,9:12:end);
tconst_mar = tconst(:,:,3:12:end);
tot_sep = tot(:,:,9:12:end);
tot_mar = tot(:,:,3:12:end);

var_sconst = var(sconst_mar,0,3,'omitnan');
var_tconst = var(tconst_mar,0,3,'omitnan');
var_tot = var(tot_mar,0,3,'omitnan');

figure
levels = 1*[0:0.05:1];
ax1 = subaxis(1,2,1,'SH',0.02,'SV',0.05);
contourf(xi,yi,(var_sconst./var_tot)',levels,'edgecolor','none')
colorbar
caxis([levels(1) levels(end)])
cmocean('deep',length(levels)-1)
ylim([-50 50])
title('% surface density variance from T, March','FontSize',6)
ax1.FontSize = 6;

ax2 = subaxis(1,2,2);
contourf(xi,yi,(var_tconst./var_tot)',levels,'edgecolor','none')
colorbar
caxis([levels(1) levels(end)])
cmocean('deep',length(levels)-1)
ylim([-50 50])
title('% surface density variance from S, March','FontSize',6)
ax2.FontSize = 6;



%% calculate wind stress curl, add to P gridded file
expname = 'interannual_tflux';
%rootdir = '/vast'; %on taris or ilko; or /batou for iter129_bulkformula, nointerannual, and interannual_southpac
rootdir = '/vast/proj/ecco'; %vast path is different on atalanta
%rootdir = '/batou';
%addpath /atalanta/home/swijffels/toolbox/seawater %on taris or ilko
addpath /home/swijffels/toolbox/seawater %on atalanta
e = matfile([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/mat/P_gridded/gridonP_ecco4r4_',expname,'.mat']);

txi = e.txi;
tyi = e.tyi;
xi = e.xi;
yi = e.yi;
yrgrid = e.yrgrid;

Re = 6378137; %radius of Earth [m]
curl_tau = NaN(length(xi),length(yi),length(yrgrid));
for tt = 1:length(yrgrid)
    [~,DtxDy] = gradient(squeeze(txi(:,:,tt)),1,1);
    [DtyDx,~] = gradient(squeeze(tyi(:,:,tt)),1,1);
    
    dr = (pi/180)*Re; %length of 1 degree lat in m
    dlat = yi(2) - yi(1);
    dlon = xi(2) - xi(1);
    
    DtxDy = DtxDy/(dr*dlat);
    Dx = dlon*dr*cosd(yi);
    
    DtyDx = DtyDx./(repmat(Dx,1,length(xi)))';
    
    curl_tau(:,:,tt) = DtyDx - DtxDy;
end

save([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/mat/P_gridded/gridonP_ecco4r4_',expname,'.mat'],'curl_tau','-append');

