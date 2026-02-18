%% CALCULATE STD, MAD, MEDIANS OF DPR DATA
% using error estimates to screen lowess fits
% _fwa = weighted linear fit (offset and slope), fwe is associated error
% _wa = simple spatial weight
addpath /atalanta/home/swijffels/work/argo/gridNSF

a = matfile('/atalanta/home/swijffels/work/argo/gridNSF/gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
xi = a.xi;
yi = a.yi;
yrgrid = a.yrgrid;
sig1grid = a.sig1grid;

[nx,ny,ns,nt]=size(a.dpr_fwa);
[dprstdo,edprmed,dprMAD,dprstd]=deal(NaN(nx,ny,ns));

if ~exist('clim')
    clim = load('/atalanta/home/swijffels/work/argo/gridNSF/sig1grid/gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end 

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/atalanta/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end
%% calculate some stats...
c = -1/(sqrt(2)*erfcinv(3/2));
for is = 1:ns
    data = squeeze(a.dpr_fwa(:,:,is,:));
    dprstdo(:,:,is) = nanstd(data,1,3);  % std of unscreened series
    edprmed(:,:,is) = nanmedian(squeeze(a.dpr_fwe(:,:,is,:)),3);  % median of error estimate in lowess
    datmed=nanmedian(data,3);
    dattemp = NaN(nx,ny,nt);
    for it = 1:nt
        dattemp(:,:,it) = abs(data(:,:,it) - datmed);
    end
    dprMAD(:,:,is) = c*nanmedian(dattemp,3);
    is
end
%% plot maps of std on sigma1 surfaces
figure
for is = select_is %1:ns
clf
pcolor(xi,yi,squeeze(dprstdo(:,:,is))'),shading flat, caxis([-1,1]*10)   % shows where there are outliers
is
pause
end

%% figure out which dpr data to screen

ibad = true(nx,ny,ns,nt);
for is = 1:ns
    mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,is)) + squeeze(clim.pr.c(1:4:end,21:4:621,is,1));
    mean_dpr = squeeze(clim.dpr.m(1:4:end,21:4:621,is)) + squeeze(clim.dpr.c(1:4:end,21:4:621,is,1));
    winter_outcrop = mean_pr <= clim_mldminmax.mldmax(1:4:end,21:4:621) | sig1grid(is) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621);

    for it = find(yrgrid>2004) %for times after beginning of 2004
        above_surface = mean_pr + squeeze(a.pr_fwa(:,:,is,it)) <= 0;
        high_error = squeeze(a.dpr_fwe(:,:,is,it)) > 3*edprmed(:,:,is); %lower from 4 to 3
        value_spike = abs(squeeze(a.dpr_fwa(:,:,is,it))) > 5*dprMAD(:,:,is); %raise this from 4 to 5
        missing_error = isnan(squeeze(a.dpr_fwe(:,:,is,it)));
        negative_thickness = mean_dpr + a.dpr_fwa(:,:,is,it) < 0;
        ibad(:,:,is,it) = winter_outcrop | above_surface | high_error | value_spike | missing_error | negative_thickness;
        it
    end
    is
end

%% screen dpr data
dpr_fwa = NaN(nx,ny,ns,nt); %this will be the final data
dpr_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns %run thru sig1 levels
    zz=squeeze(a.dpr_fwa(:,:,is,:)); %raw dpr data - nx ny nt
    zz(squeeze(ibad(:,:,is,:)))= NaN; %set dpr = NaN where ibad is true (high error, etc)
    dpr_fwa_nointerp(:,:,is,:) = zz;
    
    for ix = 1:nx % for each timeseries on this sigma1 level, linearly interpolate
        for iy = 1:ny
            %skip if not enough data points ... hopefully just over
            %land/etc
            if sum(~isnan(squeeze(zz(ix,iy,:)))) < 0.5*length(yrgrid)
                disp('too little data')
                dpr_fwa_nointerp(ix,iy,:) = NaN; %set this timeseries to NaN so it doesn't seem like "extra" non-NaNs
                continue
            end
            numbs = 1:length(yrgrid); % just a list of integers length of yrgrid
            clean_numbs = numbs(squeeze(~ibad(ix,iy,is,:))); % get rid of integers corresponding to bad data
            timeseries_dpr_fwa = squeeze(zz(ix,iy,:)); %just pulling out timeseries of dpr at this location
            clean_dpr_fwa = timeseries_dpr_fwa(~isnan(timeseries_dpr_fwa)); %get rid of NaNs
            interp_dpr_fwa = interp1(clean_numbs,clean_dpr_fwa,numbs,'linear'); %interpolate over bad points
            interp_dpr_fwa(isnan(interp_dpr_fwa)) = 0.0; %set all gaps to 0.0 instead of NaN
            interp_dpr_fwa(1:12) = NaN; %don't want anything pre-2004
            dpr_fwa(ix,iy,is,:) = interp_dpr_fwa;
        end
    end
    
    dprstd(:,:,is) = nanstd(zz,1,3); %calculate standard deviation of filtered dpr
    is

end

%% calculate fraction interpolated
%this is so dumb
dpr_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns
    zz = squeeze(a.dpr_fwa(:,:,is,:));
    for ix = 1:nx
        zzx = squeeze(zz(ix,:,:));
        for iy = 1:ny
            if sum(~isnan(squeeze(dpr_fwa(ix,iy,is,:)))) > 1 %only bother with timeseries that have been subject to interpolation
                zzxy = squeeze(zzx(iy,:));
                zzxy(squeeze(ibad(ix,iy,is,:))) = NaN;
                dpr_fwa_nointerp(ix,iy,is,:) = zzxy;
            else
                continue
            end
        end
        disp(ix)
    end
    disp(is)
end

fraction_interpolated = (sum(~isnan(dpr_fwa),'all') - sum(~isnan(dpr_fwa_nointerp),'all'))/sum(~isnan(dpr_fwa),'all');

%% plot new std (of despiked data)
for is = 1:ns %select_is %1:ns
clf
subplot(1,2,1) %download subaxis toolbox
pcolor(xi,yi,squeeze(dprstd(:,:,is))'),shading flat, caxis([-1,1]*15)
subplot(1,2,2)
pcolor(xi,yi,squeeze(dprstdo(:,:,is))'),shading flat, caxis([-1,1]*15)
is
pause
end

%% compare timeseries of dpr

ix = find(xi==100);
is = 32;

%plots different versions of dpr, along with error
for iy = 26:2:126 % within -50 to 50 latitude
    zp = squeeze(a.dpr_fwa(ix,iy,is,:)); %dpr data
    ze = squeeze(a.dpr_fwe(ix,iy,is,:)); %dpr error
    clf
    hold on
    plot(yrgrid,0*yrgrid,'k-') %just a horizontal line at zero
    plot(yrgrid,3*squeeze(edprmed(ix,iy,is))+0*yrgrid,'r-'); %plot error threshold
    plot(yrgrid,5*squeeze(dprMAD(ix,iy,is))+0*yrgrid,'b-'); %plot dpr value threshold
    plot(yrgrid,-5*squeeze(dprMAD(ix,iy,is))+0*yrgrid,'b-'); %plot dpr value threshold
    plot(yrgrid,ze,'r--') %dpr error
    h=plot(yrgrid,zp,'b',yrgrid,squeeze(a.dpr_wa(ix,iy,is,:)),'g') %blue: dpr, green: dpr standard weighted average
    set(h,'linewidth',1.5)
    ib = ibad(ix,iy,is,:);  % screen choice - here using 5 x median of the error, etc
    plot(yrgrid(ib),zp(ib),'r*') % plotting dpr where error is high
    zpnointerp = zp;
    zpnointerp(ib) = NaN;
    h=plot(yrgrid,zpnointerp,'k-') %plotting dpr where error is low
    set(h,'linewidth',2)
    zp_new = squeeze(dpr_fwa(ix,iy,is,:));
    h = plot(yrgrid,zp_new,'m-'); %plot the new, despiked, interpolated data
    legend('.','error threshold','dpr threshold','dpr threshold','dpr error','dpr_fwa','dpr_wa','dpr where ibad','dpr w low error')
    title(['latitude ',num2str(yi(iy))])
    pause
end

%% CALCULATE STD, MAD, MEDIANS OF SALINITY DATA
% _fwa = weighted linear fit (offset and slope), fwe is associated error
% _wa = simple spatial weight

addpath /home/swijffels/work/argo/gridNSF

a = matfile('/home/swijffels/work/argo/gridNSF/gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
xi = a.xi;
yi = a.yi;
yrgrid = a.yrgrid;
sig1grid = a.sig1grid;

[nx,ny,ns,nt]=size(a.sa_fwa);
[sastdo,esamed,saMAD,sastd]=deal(NaN(nx,ny,ns));

if ~exist('clim')
    clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end 

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/atalanta/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end

%% calculate some stats...
c = -1/(sqrt(2)*erfcinv(3/2));
for is = 1:ns
    data = squeeze(a.sa_fwa(:,:,is,:));
    sastdo(:,:,is) = nanstd(data,1,3);  % std of unscreened series
    esamed(:,:,is) = nanmedian(squeeze(a.sa_fwe(:,:,is,:)),3);  % median of error estimate in lowess
    datmed=nanmedian(data,3);
    dattemp = NaN(nx,ny,nt);
    for it = 1:nt
        dattemp(:,:,it) = abs(data(:,:,it) - datmed);
    end
    saMAD(:,:,is) = c*nanmedian(dattemp,3);  
    is
end
%% figure out which salinity data to screen

ibad = true(nx,ny,ns,nt);
for is = 1:ns
    mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,is)) + squeeze(clim.pr.c(1:4:end,21:4:621,is,1));
    mean_sa = squeeze(clim.sa.m(1:4:end,21:4:621,is)) + squeeze(clim.sa.c(1:4:end,21:4:621,is,1));
    winter_outcrop = mean_pr <= clim_mldminmax.mldmax(1:4:end,21:4:621) | sig1grid(is) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621);

    for it = find(yrgrid>2004) %for times after beginning of 2004
        above_surface = mean_pr + squeeze(a.pr_fwa(:,:,is,it)) <= 0;
        high_error = squeeze(a.sa_fwe(:,:,is,it)) > 3*esamed(:,:,is);
        value_spike = abs(squeeze(a.sa_fwa(:,:,is,it))) > 5*saMAD(:,:,is);
        missing_error = isnan(squeeze(a.sa_fwe(:,:,is,it)));
        negative_salinity = mean_sa + a.sa_fwa(:,:,is,it) < 0;
        ibad(:,:,is,it) = winter_outcrop | above_surface | high_error | missing_error | negative_salinity;
    it
    end
end

%% screen salinity data
sa_fwa = NaN(nx,ny,ns,nt); %this will be the final data
sa_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns %run thru sig1 levels
    zz=squeeze(a.sa_fwa(:,:,is,:)); %raw sa data - nx ny nt
    zz(squeeze(ibad(:,:,is,:)))= NaN; %set sa = NaN where ibad is true (high error, etc)
    sa_fwa_nointerp(:,:,is,:) = zz;
    
    for ix = 1:nx % for each timeseries on this sigma1 level, linearly interpolate
        for iy = 1:ny
            %skip if not enough data points ... hopefully just over
            %land/etc
            if sum(~isnan(squeeze(zz(ix,iy,:)))) < 0.5*length(yrgrid)
                disp('too little data')
                sa_fwa_nointerp(ix,iy,:) = NaN; %set this timeseries to NaN so it doesn't seem like "extra" non-NaNs
                continue
            end
            numbs = 1:length(yrgrid); % just a list of integers length of yrgrid
            clean_numbs = numbs(squeeze(~ibad(ix,iy,is,:))); % get rid of integers corresponding to bad data
            timeseries_sa_fwa = squeeze(zz(ix,iy,:)); %just pulling out timeseries of sa at this location
            clean_sa_fwa = timeseries_sa_fwa(~isnan(timeseries_sa_fwa)); %get rid of NaNs
            interp_sa_fwa = interp1(clean_numbs,clean_sa_fwa,numbs,'linear'); %interpolate over bad points
            interp_sa_fwa(isnan(interp_sa_fwa)) = 0.0; %set all gaps to 0.0 instead of NaN
            interp_sa_fwa(1:12) = NaN; %don't want anything pre-2004
            sa_fwa(ix,iy,is,:) = interp_sa_fwa;
        end
    end
    
    sastd(:,:,is) = nanstd(zz,1,3); %calculate standard deviation of filtered sa
    is

end

%% calculate fraction interpolated
%this is so dumb
sa_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns
    zz = squeeze(a.sa_fwa(:,:,is,:));
    for ix = 1:nx
        zzx = squeeze(zz(ix,:,:));
        for iy = 1:ny
            if sum(~isnan(squeeze(sa_fwa(ix,iy,is,:)))) > 1 %only bother with timeseries that have been subject to interpolation
                zzxy = squeeze(zzx(iy,:));
                zzxy(squeeze(ibad(ix,iy,is,:))) = NaN;
                sa_fwa_nointerp(ix,iy,is,:) = zzxy;
            else
                continue
            end
        end
        disp(ix)
    end
    disp(is)
end

fraction_interpolated = (sum(~isnan(sa_fwa),'all') - sum(~isnan(sa_fwa_nointerp),'all'))/sum(~isnan(sa_fwa),'all');

%% compare timeseries of salinity

ix = find(xi==200);
is = 38;

%plots different versions of salinity, along with error
for iy = 26:2:126 % within -50 to 50 latitude
    zp = squeeze(a.sa_fwa(ix,iy,is,:)); %sa data
    ze = squeeze(a.sa_fwe(ix,iy,is,:)); %sa error
    clf
    hold on
    plot(yrgrid,0*yrgrid,'k-') %just a horizontal line at zero
    plot(yrgrid,3*squeeze(esamed(ix,iy,is))+0*yrgrid,'r-'); %plot error threshold
    plot(yrgrid,4*squeeze(saMAD(ix,iy,is))+0*yrgrid,'b-'); %plot salinity value threshold
    plot(yrgrid,-4*squeeze(saMAD(ix,iy,is))+0*yrgrid,'b-'); %plot salinity value threshold
    plot(yrgrid,ze,'r--') %salinity error
    h=plot(yrgrid,zp,'b',yrgrid,squeeze(a.sa_wa(ix,iy,is,:)),'g') %blue: salinity, green: salinity standard weighted average
    set(h,'linewidth',1.5)
    ib = ibad(ix,iy,is,:);  % screen choice - here using 5 x median of the error, etc
    plot(yrgrid(ib),zp(ib),'r*') % plotting salinity where error is high
    zpnointerp = zp;
    zpnointerp(ib) = NaN;
    h=plot(yrgrid,zpnointerp,'k-') %plotting salinity where error is low
    set(h,'linewidth',2)
    zp_new = squeeze(sa_fwa(ix,iy,is,:));
    h = plot(yrgrid,zp_new,'m-'); %plot the new, despiked, interpolated data
    legend('.','error threshold','sa threshold','sa threshold','sa error','sa_fwa','sa_wa','sa where ibad','sa w low error','location','southoutside')
    title(['latitude ',num2str(yi(iy))])
    pause
end

%% CALCULATE STD, MAD, MEDIANS OF PRESSURE DATA
% _fwa = weighted linear fit (offset and slope), fwe is associated error
% _wa = simple spatial weight

addpath /atalanta/home/swijffels/work/argo/gridNSF

a = matfile('/atalanta/home/swijffels/work/argo/gridNSF/gridonSigma1_anomalies_Argo_CORA_sla_superobs_Sep2022_1998_2020.mat');
xi = a.xi;
yi = a.yi;
yrgrid = a.yrgrid;
sig1grid = a.sig1grid;

[nx,ny,ns,nt]=size(a.pr_fwa);
[prstdo,eprmed,prMAD,prstd]=deal(NaN(nx,ny,ns));

if ~exist('clim')
    clim = load('gridonSig1_climatology_cora_argo_huber25_sla_ug_lscovtar_200.mat');
end 

if ~exist('clim_mldminmax')
    clim_mldminmax = load('/atalanta/home/swijffels/work/argo/gridNSF/climfiles/clim_mldminmax.mat');
end

%% calculate some stats...
c = -1/(sqrt(2)*erfcinv(3/2));
for is = 1:ns
    data = squeeze(a.pr_fwa(:,:,is,:));
    prstdo(:,:,is) = nanstd(data,1,3);  % std of unscreened series
    eprmed(:,:,is) = nanmedian(squeeze(a.pr_fwe(:,:,is,:)),3);  % median of error estimate in lowess
    datmed=nanmedian(data,3);
    dattemp = NaN(nx,ny,nt);
    for it = 1:nt
        dattemp(:,:,it) = abs(data(:,:,it) - datmed);
    end
    prMAD(:,:,is) = c*nanmedian(dattemp,3);  
    is
end
%% figure out which pressure data to screen

ibad = true(nx,ny,ns,nt);
for is = 1:ns
    mean_pr = squeeze(clim.pr.m(1:4:end,21:4:621,is)) + squeeze(clim.pr.c(1:4:end,21:4:621,is,1));
    winter_outcrop = mean_pr <= clim_mldminmax.mldmax(1:4:end,21:4:621) | sig1grid(is) <= clim_mldminmax.mlsig1max(1:4:end,21:4:621);

    for it = find(yrgrid>2004) %for times after beginning of 2004
        above_surface = mean_pr + squeeze(a.pr_fwa(:,:,is,it)) <= 0;
        high_error = squeeze(a.pr_fwe(:,:,is,it)) > 3*eprmed(:,:,is);
        value_spike = abs(squeeze(a.pr_fwa(:,:,is,it))) > 5*prMAD(:,:,is);
        missing_error = isnan(squeeze(a.pr_fwe(:,:,is,it)));
        ibad(:,:,is,it) = winter_outcrop | above_surface | high_error | missing_error;
    it
    end
end

%% screen pressure data
pr_fwa = NaN(nx,ny,ns,nt); %this will be the final data
pr_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns %run thru sig1 levels
    zz=squeeze(a.pr_fwa(:,:,is,:)); %raw pr data - nx ny nt
    zz(squeeze(ibad(:,:,is,:)))= NaN; %set pr = NaN where ibad is true (high error, etc)
    pr_fwa_nointerp(:,:,is,:) = zz;
    
    for ix = 1:nx % for each timeseries on this sigma1 level, linearly interpolate
        for iy = 1:ny
            %skip if not enough data points ... hopefully just over
            %land/etc
            if sum(~isnan(squeeze(zz(ix,iy,:)))) < 0.5*length(yrgrid)
                disp('too little data')
                pr_fwa_nointerp(ix,iy,:) = NaN; %set this timeseries to NaN so it doesn't seem like "extra" non-NaNs
                continue
            end
            numbs = 1:length(yrgrid); % just a list of integers length of yrgrid
            clean_numbs = numbs(squeeze(~ibad(ix,iy,is,:))); % get rid of integers corresponding to bad data
            timeseries_pr_fwa = squeeze(zz(ix,iy,:)); %just pulling out timeseries of pr at this location
            clean_pr_fwa = timeseries_pr_fwa(~isnan(timeseries_pr_fwa)); %get rid of NaNs
            interp_pr_fwa = interp1(clean_numbs,clean_pr_fwa,numbs,'linear'); %interpolate over bad points
            interp_pr_fwa(isnan(interp_pr_fwa)) = 0.0; %set all gaps to 0.0 instead of NaN
            interp_pr_fwa(1:12) = NaN; %don't want anything pre-2004
            pr_fwa(ix,iy,is,:) = interp_pr_fwa;
        end
    end
    
    prstd(:,:,is) = nanstd(zz,1,3); %calculate standard deviation of filtered pr
    is

end

%% calculate fraction interpolated
%this is so dumb
pr_fwa_nointerp = NaN(nx,ny,ns,nt);

for is = 1:ns
    zz = squeeze(a.pr_fwa(:,:,is,:));
    for ix = 1:nx
        zzx = squeeze(zz(ix,:,:));
        for iy = 1:ny
            if sum(~isnan(squeeze(pr_fwa(ix,iy,is,:)))) > 1 %only bother with timeseries that have been subject to interpolation
                zzxy = squeeze(zzx(iy,:));
                zzxy(squeeze(ibad(ix,iy,is,:))) = NaN;
                pr_fwa_nointerp(ix,iy,is,:) = zzxy;
            else
                continue
            end
        end
        disp(ix)
    end
    disp(is)
end

fraction_interpolated = (sum(~isnan(pr_fwa),'all') - sum(~isnan(pr_fwa_nointerp),'all'))/sum(~isnan(pr_fwa),'all');

%% compare timeseries of pressure

ix = find(xi==270);
is = 38;

%plots different versions of pressure, along with error
for iy = 26:2:126 % within -50 to 50 latitude
    zp = squeeze(a.pr_fwa(ix,iy,is,:)); %pr data
    ze = squeeze(a.pr_fwe(ix,iy,is,:)); %pr error
    clf
    hold on
    plot(yrgrid,0*yrgrid,'k-') %just a horizontal line at zero
    plot(yrgrid,3*squeeze(eprmed(ix,iy,is))+0*yrgrid,'r-'); %plot error threshold
    plot(yrgrid,5*squeeze(prMAD(ix,iy,is))+0*yrgrid,'b-'); %plot pressure value threshold
    plot(yrgrid,-5*squeeze(prMAD(ix,iy,is))+0*yrgrid,'b-'); %plot pressure value threshold
    plot(yrgrid,ze,'r--') %pressure error
    h=plot(yrgrid,zp,'b',yrgrid,squeeze(a.pr_wa(ix,iy,is,:)),'g') %blue: pressure, green: pressure standard weighted average
    set(h,'linewidth',1.5)
    ib = ibad(ix,iy,is,:);  % screen choice - here using 5 x median of the error, etc
    plot(yrgrid(ib),zp(ib),'r*') % plotting pressure where error is high
    zpnointerp = zp;
    zpnointerp(ib) = NaN;
    h=plot(yrgrid,zpnointerp,'k-') %plotting pressure where error is low
    set(h,'linewidth',2)
    zp_new = squeeze(pr_fwa(ix,iy,is,:));
    h = plot(yrgrid,zp_new,'m-'); %plot the new, despiked, interpolated data
    legend('.','error threshold','pr threshold','pr threshold','pr error','pr_fwa','pr_wa','pr where ibad','pr w low error','location','southoutside')
    title(['latitude ',num2str(yi(iy))])
    pause
end
