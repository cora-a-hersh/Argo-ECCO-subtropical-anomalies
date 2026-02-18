%% takes postprocessed ECCO output and converts to .mat files for pressure- and sig1-gridded data
%what is the name of the ECCO experiment whose output you want to put into
%.mat files?

expname = 'interannual_eqind';

%on taris or ilko
%rootdir = '/vast';
%rootdir = '/batou'; %for interannual_southpac, nointerannual, and
%iter129_bulkformula
%on atalanta
rootdir = '/vast/proj/ecco';
%rootdir = '/batou';

%are vector fields (e.g. velocity) available? 1 = yes, 0 = no
if exist([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/NVELMASS'],'dir')
    vel_avail = 1;
else
    vel_avail = 0;
end

addpath /home/swijffels/toolbox/seawater
addpath /home/swijffels/toolbox/csirolib
addpath /home/swijffels/work/argo/matlab

%{
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
%}

%after running GCM, two steps: regularpoles, sigma1 gridding
%scripts directory on github

datdir = [rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles'];
savedatdir = [rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/mat/'];

%tname = [datdir,'/THETA/THETA_on_sigma1_1996_02.nc'];
nyears = 2017-1992 +1;  
nt = nyears*12;  % year range

%% keep sigma1-gridded files on native grid:
tname = [datdir,'/THETA/THETA_on_sigma1_1996_02.nc'];
        x = ncread(tname,'lon');
        x = x(:);
           
        xs = rem(x+360,360);
        [xs,isx]=sort(xs);
        
%         ibad = abs(xs - 180) < .5;
%         isx(ibad) = [];  % bad data interpolation near dateline
%         xs(ibad) = [];
%        xs = [-.25;xs;360.25];
%        isx = [find(x==-.25);isx;find(x==.25)];
        
        y = ncread(tname,'lat');
        y  = y(:);
        sig1grid = ncread(tname,'sigma-1'); % kgm-3, name changed to sigma-1 from sigma1 in Nov 2021
        xi=xs;yi = y;
 
 nx = length(xi);ny=length(yi); nz=length(sig1grid);
 
%% just add DRHODR (drho/dz)

for iyear = 1993:2017 %1993:2017 
    tic
    cyear = int2str(iyear)
    [drhodri]=deal(NaN(nx,ny,nz,12));   % cannot do all years  in one go! Not enough memory on atalanta: do in 12 months blocks
    
    dnum = [];
    
    for it = 1:12
        cmon = int2str(it+100);
        cmon(1) = [];
     drhodrname = [datdir,'/DRHODR/DRHODR_on_sigma1_',cyear,'_',cmon,'.nc']; 
     mday = datenum(iyear,it,15,12,0,0); % middle of the month
     dnum = [dnum;mday];  
     drhodr1 = ncread(drhodrname,'DRHODR');
        
         for iz = 1:nz
            drhodri(:,:,iz,it) = drhodr1(isx,:,iz);
         end
     
       datestr(mday)   

    end
 
  yrgrid = iyear + [0.5:1:(12)]/12;

savename = [savedatdir,'sig1_gridded/','drhodr_gridonSig1_ecco4r4_',expname,'_',cyear,'.mat']

save -v7.3 junk.mat drhodri

eval(['!mv junk.mat ',savename]);
toc

end
%% concatenate yearly drhodr files
tic
addpath([savedatdir,'sig1_gridded']);
yrfiles = dir([savedatdir,'sig1_gridded/drhodr*']);
nyr = length(yrfiles); %number of years of data
yrnames = cell(size(yrfiles));
for ii = 1:length(yrfiles) %probably a better way to do this but w/e
    yrnames{ii} = yrfiles(ii).name;
end

load(yrnames{3})

%initialize larger total-data arrays
drhodri_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));

for yy = 1:nyr %for each year file
   load(yrnames{yy}) %load in that year of data
   drhodri_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = drhodri;
   clear drhodri
end

drhodri = drhodri_tot;

%save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'drhodri','-append')
save(['/vast/proj/ecco/ECCOv4r4/exps/interannual_northpac/run/regularpoles/mat1/sig1_gridded/drhodr_gridonSig1_ecco4r4_',expname,'.mat'],'drhodri')

toc
  
%% convert sigma1-gridded files to .mat files

if ~exist([savedatdir,'sig1_gridded'],'dir')
    mkdir([savedatdir,'sig1_gridded'])
end

for iyear = 1993:2017
    tic
    cyear = int2str(iyear)
    [si,pri,pti,vi,ui,drhodri]=deal(NaN(nx,ny,nz,12));   % cannot do all years  in one go! Not enough memory on atalanta: do in 12 months blocks
    
    dnum = [];
    
    for it = 1:12
        cmon = int2str(it+100);
        cmon(1) = [];
     tname = [datdir,'/THETA/THETA_on_sigma1_',cyear,'_',cmon,'.nc'];
     sname = [datdir,'/SALT/SALT_on_sigma1_',cyear,'_',cmon,'.nc']; 
     if vel_avail == 1
        uname = [datdir,'/EVELMASS/EVELMASS_on_sigma1_',cyear,'_',cmon,'.nc']; 
        vname = [datdir,'/NVELMASS/NVELMASS_on_sigma1_',cyear,'_',cmon,'.nc']; 
     end
     pname = [datdir,'/p/p_on_sigma1_',cyear,'_',cmon,'.nc'];  
     drhodrname = [datdir,'/DRHODR/DRHODR_on_sigma1_',cyear,'_',cmon,'.nc']; 
     mday = datenum(iyear,it,15,12,0,0); % middle of the month
     dnum = [dnum;mday];  
     
     pt1 = ncread(tname,'THETA') ;
     s1 = ncread(sname,'SALT') ; 
     if vel_avail == 1
        u1 = ncread(uname,'EVELMASS') ;  % m/s
        v1 = ncread(vname,'NVELMASS') ; 
     end
     pr1 = ncread(pname,'p') ; % dbar
     drhodr1 = ncread(drhodrname,'DRHODR'); %kg/m^4
        
         for iz = 1:nz
%           junk = squeeze( pr1(:,:,iz));  % depth is offset by one grid
%           point and not yet fixed. 
%           ibad = junk < 0 | junk > di;
%           bmask = ones(nx,ny);
%           bmask(ibad) = NaN;
%           pri(:,:,iz,it) =  pr1(isx,:,iz).*bmask(isx,:);
%           si(:,:,iz,it) =  s1(isx,:,iz).*bmask(isx,:);
%           pti(:,:,iz,it) =  pt1(isx,:,iz).*bmask(isx,:);
%           ui(:,:,iz,it) =  u1(isx,:,iz).*bmask(isx,:);
%           vi(:,:,iz,it) =  v1(isx,:,iz).*bmask(isx,:);

            pri(:,:,iz,it) =  pr1(isx,:,iz);
            si(:,:,iz,it) =  s1(isx,:,iz);
            pti(:,:,iz,it) =  pt1(isx,:,iz);
            drhodri(:,:,iz,it) = drhodr1(isx,:,iz);
            if vel_avail == 1
                ui(:,:,iz,it) =  u1(isx,:,iz);
                vi(:,:,iz,it) =  v1(isx,:,iz);
            end
         end
     
       datestr(mday)   

    end
 
  yrgrid = iyear + [0.5:1:(12)]/12;

readme = 'made from J Gebbie & C Hersh run of MITgcm for eccov4r4 by convert_ECCO_nc_to_mat.m  Jan 2025';

savename = [savedatdir,'sig1_gridded/','gridonSig1_ecco4r4_',expname,'_',cyear,'.mat']
if vel_avail == 1
    save -v7.3 junk.mat xi yi yrgrid si pti pri drhodri ui vi sig1grid dnum yrgrid readme 
else
    save -v7.3 junk.mat xi yi yrgrid si pti pri drhodri sig1grid dnum yrgrid readme 
end
eval(['!mv junk.mat ',savename]);
toc

end

%% now concatenate all separate year files into one larger file
%go variable by variable to see if it stops this from crashing...

addpath([savedatdir,'sig1_gridded']);
yrfiles = dir([savedatdir,'sig1_gridded']);
yrfiles = yrfiles(3:end); %cut off . and ..
nyr = length(yrfiles); %number of years of data
yrnames = cell(size(yrfiles));
for ii = 1:length(yrfiles) %probably a better way to do this but w/e
    yrnames{ii} = yrfiles(ii).name;
end

load(yrnames{3})

varnames = {'pri','pti','si','ui','vi','drhodri','dnum','yrgrid'};

for vv = 1:length(varnames) %for each variable to concatenate
    disp(vv)
    if vv == 1
        dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
        for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = pri;
        end
        pri = dat_tot;
        save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'pri','sig1grid','xi','yi','-v7.3')
        clear pri dat_tot
    elseif vv == 2
       dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = pti;
       end
       pti = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'pti','-append')
       clear pti dat_tot
    elseif vv == 3
       dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = si;
       end
       si = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'si','-append')
       clear si dat_tot
    elseif vv == 4
       dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = ui;
       end
       ui = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'ui','-append')
       clear ui dat_tot
    elseif vv == 5
       dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = vi;
       end
       vi = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'vi','-append')
       clear vi dat_tot
    elseif vv == 6
       dat_tot = NaN(length(xi),length(yi),length(sig1grid),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = drhodri;
       end
       drhodri = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'drhodri','-append')
       clear drhodri dat_tot
    elseif vv == 7
       dat_tot = NaN(nyr*length(yrgrid),1);
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot((1+(yy-1)*12):(yy*12)) = dnum;
       end
       dnum = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'dnum','-append')
       clear dnum dat_tot
    elseif vv == 8
       dat_tot = NaN(1,nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot((1+(yy-1)*12):(yy*12)) = yrgrid;
       end
       yrgrid = dat_tot;
       save([savedatdir,'sig1_gridded/gridonSig1_ecco4r4_',expname,'.mat'],'yrgrid','-append')
       clear yrgrid dat_tot
    end
end
%take out overall mean
%calculate mean from 2004-2017 to be more comparable to Argo
%could go back into old code and calculate mean from 2004-2017 as well

%% keep pressure-gridded files on native grid:

tname = [datdir,'/SALT/SALT_1992_01.nc'];
x = ncread(tname,'lon');
x = x(:);

xs = rem(x+360,360);
[xs,isx]=sort(xs);

y = ncread(tname,'lat');
y = y(:);
z = ncread(tname,'depth');
z = z(:);
xi=xs; yi=y; zi=z;
 
nx=length(xi);ny=length(yi); nz=length(zi);

%% convert pressure-gridded files to .mat files
addpath([rootdir,'/ECCOv4r4/exps/',expname,'/run/regularpoles/'])

if ~exist([savedatdir,'P_gridded'],'dir')
    mkdir([savedatdir,'P_gridded'])
end

for iyear = 1993:2017
    tic
    cyear = int2str(iyear)
    [si,pdi,pti,vi,ui,drhodri]=deal(NaN(nx,ny,nz,12)); % cannot do all years in one go! Not enough memory on atalanta: do in 12 months blocks
    [mxi,fwi,qi,txi,tyi,tflxi,qswi] = deal(NaN(nx,ny,12));
    %pdi: potential density ref to surface
    %mxi: mixed layer depth - 2D
    %fwi: fresh water flux - 2D
    %qi: Qnet surface heat flux - 2D
    %tflxi: TFLUX - 2D
    %qswi: short wave heat flux - 2D
    %txi: tau x, wind - 2D
    %tyi: tau y, wind - 2D
    dnum = [];
    
    for it = 1:12
        cmon = int2str(it+100);
        cmon(1) = [];
        sname = [datdir,'/SALT/SALT_',cyear,'_',cmon,'.nc']; 
        dname = [datdir,'/RHOAnoma/RHOAnoma_',cyear,'_',cmon,'.nc']; %potential density
        tname = [datdir,'/THETA/THETA_',cyear,'_',cmon,'.nc'];
        if vel_avail == 1
            uname = [datdir,'/EVELMASS/EVELMASS_',cyear,'_',cmon,'.nc']; 
            vname = [datdir,'/NVELMASS/NVELMASS_',cyear,'_',cmon,'.nc']; 
        end
        drhname = [datdir,'/DRHODR/DRHODR_',cyear,'_',cmon,'.nc'];
        mxname = [datdir,'/MXLDEPTH/MXLDEPTH_',cyear,'_',cmon,'.nc']; 
        fwname = [datdir,'/oceFWflx/oceFWflx_',cyear,'_',cmon,'.nc']; 
        qname = [datdir,'/oceQnet/oceQnet_',cyear,'_',cmon,'.nc']; 
        txname = [datdir,'/oceTAUE/oceTAUE_',cyear,'_',cmon,'.nc']; 
        tyname = [datdir,'/oceTAUN/oceTAUN_',cyear,'_',cmon,'.nc']; 
        tflxname = [datdir,'/TFLUX/TFLUX_',cyear,'_',cmon,'.nc']; 
        qswname = [datdir,'/oceQsw/oceQsw_',cyear,'_',cmon,'.nc']; 
        mday = datenum(iyear,it,15,12,0,0); % middle of the month
        dnum = [dnum;mday];  
     
        s1 = ncread(sname,'SALT'); 
        d1 = ncread(dname,'RHOAnoma');
        pt1 = ncread(tname,'THETA');
        u1 = ncread(uname,'EVELMASS');  % m/s
        v1 = ncread(vname,'NVELMASS'); 
        drh1 = ncread(drhname,'DRHODR');
        mx1 = ncread(mxname,'MXLDEPTH');
        fw1 = ncread(fwname,'oceFWflx');
        q1 = ncread(qname,'oceQnet');
        tx1 = ncread(txname,'oceTAUE');
        ty1 = ncread(tyname,'oceTAUN');
        tflx1 = ncread(tflxname,'TFLUX');
        qsw1 = ncread(qswname,'oceQsw');
        
        for iz = 1:nz
          si(:,:,iz,it) =  s1(isx,:,iz);
          pdi(:,:,iz,it) = d1(isx,:,iz);
          pti(:,:,iz,it) =  pt1(isx,:,iz);
          ui(:,:,iz,it) =  u1(isx,:,iz);
          vi(:,:,iz,it) =  v1(isx,:,iz);
          drhodri(:,:,iz,it) = drh1(isx,:,iz);
          mxi(:,:,it) = mx1(isx,:);
          fwi(:,:,it) = fw1(isx,:);
          qi(:,:,it) = q1(isx,:);
          txi(:,:,it) = tx1(isx,:);
          tyi(:,:,it) = ty1(isx,:);
          tflxi(:,:,it) = tflx1(isx,:);
          qswi(:,:,it) = qsw1(isx,:);
          
         end
     
       datestr(mday)   

    end
    
  yrgrid = iyear + [0.5:1:(12)]/12;

readme = 'made from J Gebbie and C Hersh run of MITgcm for eccov4r4 by convert_ECCO_nc_to_mat.m Dec 2024';

savename = [savedatdir,'P_gridded/','gridonP_ecco4r4_',expname,'_',cyear,'.mat']

save -v7.3 junk.mat xi yi zi yrgrid si pdi pti mxi fwi qi txi tyi tflxi qswi ui vi drhodri dnum yrgrid readme

eval(['!mv junk.mat ',savename]); 
toc

end

%% now concatenate all separate year files into one larger file

addpath([savedatdir,'P_gridded']);
yrfiles = dir([savedatdir,'P_gridded']);
yrfiles = yrfiles(3:end); %cut off . and ..
nyr = length(yrfiles); %number of years of data
yrnames = cell(size(yrfiles));
for ii = 1:length(yrfiles) %probably a better way to do this but w/e
    yrnames{ii} = yrfiles(ii).name;
end

load(yrnames{3})

varnames = {'pdi','pti','si','ui','vi','drhodri','mxi','txi','tyi','qi','fwi','tflxi','qswi','dnum','yrgrid'};

for vv = 1:length(varnames) %for each variable to concatenate
    disp(vv)
    if vv == 1
        dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
        for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = pdi;
        end
        pdi = dat_tot;
        save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'pdi','zi','xi','yi','-v7.3')
        clear pdi dat_tot
    elseif vv == 2
       dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = pti;
       end
       pti = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'pti','-append')
       clear pti dat_tot
    elseif vv == 3
       dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = si;
       end
       si = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'si','-append')
       clear si dat_tot
    elseif vv == 4
       dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = ui;
       end
       ui = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'ui','-append')
       clear ui dat_tot
    elseif vv == 5
       dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = vi;
       end
       vi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'vi','-append')
       clear vi dat_tot
    elseif vv == 6
       dat_tot = NaN(length(xi),length(yi),length(zi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,:,(1+(yy-1)*12):(yy*12)) = drhodri;
       end
       drhodri = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'drhodri','-append')
       clear drhodri dat_tot
    elseif vv == 7
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = mxi;
       end
       mxi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'mxi','-append')
       clear mxi dat_tot
    elseif vv == 8
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = txi;
       end
       txi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'txi','-append')
       clear txi dat_tot
    elseif vv == 9
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = tyi;
       end
       tyi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'tyi','-append')
       clear tyi dat_tot
    elseif vv == 10
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = qi;
       end
       qi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'qi','-append')
       clear qi dat_tot
    elseif vv == 11
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = fwi;
       end
       fwi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'fwi','-append')
       clear fwi dat_tot
    elseif vv == 12
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = tflxi;
       end
       tflxi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'tflxi','-append')
       clear tflxi dat_tot
    elseif vv == 13
       dat_tot = NaN(length(xi),length(yi),nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot(:,:,(1+(yy-1)*12):(yy*12)) = qswi;
       end
       qswi = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'qswi','-append')
       clear qswi dat_tot
    elseif vv == 14
       dat_tot = NaN(nyr*length(yrgrid),1);
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot((1+(yy-1)*12):(yy*12)) = dnum;
       end
       dnum = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'dnum','-append')
       clear dnum dat_tot
    elseif vv == 15
       dat_tot = NaN(1,nyr*length(yrgrid));
       for yy = 1:nyr %for each year file
            load(yrnames{yy},varnames{vv}); %load in that year of data
            dat_tot((1+(yy-1)*12):(yy*12)) = yrgrid;
       end
       yrgrid = dat_tot;
       save([savedatdir,'P_gridded/gridonP_ecco4r4_',expname,'.mat'],'yrgrid','-append')
       clear yrgrid dat_tot
    end
end
