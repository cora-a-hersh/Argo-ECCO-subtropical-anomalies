% Cora Hersh, updated Oct 17 2023
load /home/chersh/SpiceAnomalies/CheltonRossbyAtlas.mat;
addpath /home/chersh/SpiceAnomalies/strmat_sla_200rmax;
addpath /home/swijffels/toolbox/seawater;

%% turn Chelton Rossby atlas in to 2D arrays instead of lists of points
lat = -75.5:1:89.5;
lon = 0.5:1:359.5;

phspeed = NaN(length(lat),length(lon));
defrad = NaN(length(lat),length(lon));

for ii = 1:length(lat)
    for jj = 1:length(lon)
        idx = CheltonAtlas.latitude == lat(ii) & CheltonAtlas.longitude== lon(jj);
        if sum(idx) ~= 0
           phspeed(ii,jj) = CheltonAtlas.phase_speed(idx);
           defrad(ii,jj) = CheltonAtlas.radius_deformation(idx);
        end
    end
end

%% interpolate to 1/8 degree (to sort of match streamline location)

lonwrap = -0.5:1:360.5; %need to be able to interpolate down to 0 and up to 360 degrees
[Lonwrap, Lat] = meshgrid(lonwrap,lat);

phspeedwrap = [phspeed(:,end),phspeed(:,:),phspeed(:,1)]; % to match lonwrap
defradwrap = [defrad(:,end),defrad(:,:),defrad(:,1)]; % to match lonwrap

latfine = -75.5:0.125:89.5;
lonfine = 0:0.125:360;

[Lonfine,Latfine] = meshgrid(lonfine,latfine);

phspeedfine = interp2(Lonwrap,Lat,phspeedwrap,Lonfine,Latfine);
defradfine = interp2(Lonwrap,Lat,defradwrap,Lonfine,Latfine);

%save('/atalanta/home/chersh/SpiceAnomalies/CheltonRossbyAtlas.mat','phspeedfine','defradfine','latfine','lonfine','-append');

%% calculate rossby wave phase speeds
rossbyspeed = NaN(length(lon),length(lat));
for ii = 1:length(lat)
    for jj = 1:length(lon)
        beta = 2*7.292*10^(-5)*cosd(lat(ii))/(6371*10^3);
        f = 2*7.292*10^(-5)*sind(lat(ii));
        rossbyspeed(jj,ii) = beta*(phspeed(ii,jj)/f)^2;
    end
end

rossbyspeed(:,71:82) = phspeed(71:82,:)'./3;
levels = [0,0.01,0.02,0.04,0.06,0.08,0.10,0.15,0.2,0.3];
figure
cont = contourf(lon,lat,rossbyspeed',levels);
c1 = colorbar;
caxis([levels(1) levels(end)])
c1.Label.String = 'm/s';
title('1st mode Rossby phase speeds calculated from Chelton 1998')
%% add estimated Rossby wave phase speeds to streamline .mat files
cd /home/chersh/SpiceAnomalies/strmat_sla_200rmax;
strmfiles = dir('/home/chersh/SpiceAnomalies/strmat_sla_200rmax');
strmfiles = strmfiles(3:end); %cut off . and ..
strmnames = cell(size(strmfiles));
for ii = 1:length(strmfiles) %probably a better way to do this but w/e
    strmnames{ii} = strmfiles(ii).name;
end

for kk = 1:length(strmfiles)
    load(strmnames{kk})
    disp(kk)
    if exist('ostrm') %add rossby speeds to argo streamlines
        for ii = 1:length(ostrm) %for each streamline
            longitude = ostrm(ii).lon;
            latitude = ostrm(ii).lat;
            %deal with points with negative or 360+ longitudes
            if any(longitude > 360)
                longitude(longitude > 360) = longitude(longitude>360)-360;
            end
            if any(longitude < 0)
                longitude(longitude < 0) = longitude(longitude<0)+360;
            end
            rossbyspeed = NaN(size(ostrm(ii).speed));
            rossbyspeed_plusU = NaN(size(ostrm(ii).speed));
            for jj = 1:length(longitude(~isnan(longitude) & ~isnan(latitude))) %for each point on the streamline
                if latitude(jj) > max(latfine) || latitude(jj) < min(latfine)
                    continue %can't interpolate to higher latitudes than Chelton map
                end
                closelat = abs(latfine - latitude(jj));
                [mmlat,idxlat] = min(closelat); %index of closest latitude grid to streamline point
                closelon = abs(lonfine - longitude(jj));
                [mmlon,idxlon] = min(closelon); %index of closest longitude grid to streamline point
                if latitude(jj) < 10 && latitude(jj) > -10 %use first-mode equatorial rossby phase speed for near equator
                    rossbyspeed(jj) = phspeedfine(idxlat,idxlon)/3;
                else %outside equatorial region
                    %calculate beta and f at this latitude
                    beta = 2*7.292*10^(-5)*cosd(latitude(jj))/(6371*10^3);
                    f = 2*7.292*10^(-5)*sind(latitude(jj));
                    rossbyspeed(jj) = beta*(phspeedfine(idxlat,idxlon)/f)^2;
                end
                
                if jj == length(ostrm(ii).lat)
                    rossbyspeed(jj) = rossbyspeed(jj-1); %just make last speed same as previous...
                    rossbyspeed_plusU(jj) = rossbyspeed(jj-1) - ostrm(ii).u(jj-1);
                else
                    %add background -u velocity
                    rossbyspeed_plusU(jj) = rossbyspeed(jj) - ostrm(ii).u(jj);
                    %project speed onto west direction
                    [~,phaseangle] = sw_dist([latitude(jj),latitude(jj+1)],[longitude(jj),longitude(jj+1)]);
                    rossbyspeed(jj) = -rossbyspeed(jj)*cosd(phaseangle); %or x -1 instead of +180 degrees
                    rossbyspeed_plusU(jj) = -rossbyspeed_plusU(jj)*cosd(phaseangle);
                end
                
            end
            ostrm(ii).rossbyspeed = rossbyspeed;
            ostrm(ii).rossbyspeed_plusU = rossbyspeed_plusU;
        end
        save(strmnames{kk},'ostrm');
    elseif exist('ostrme') %do the same for the ecco streamlines
        for ii = 1:length(ostrme) %for each streamline
            longitude = ostrme(ii).lon;
            latitude = ostrme(ii).lat;
            if any(longitude > 360)
                longitude(longitude > 360) = longitude(longitude>360)-360;
            end
            if any(longitude < 0)
                longitude(longitude < 0) = longitude(longitude<0)+360;
            end
            rossbyspeed = NaN(size(ostrme(ii).speed));
            rossbyspeed_plusU = NaN(size(ostrme(ii).speed));
            for jj = 1:length(longitude(~isnan(longitude))) %for each point on the streamline
                if latitude(jj) > max(latfine) || latitude(jj) < min(latfine)
                    continue
                end
                closelat = abs(latfine - latitude(jj));
                [mmlat,idxlat] = min(closelat); %index of closest latitude grid to streamline point
                closelon = abs(lonfine - longitude(jj));
                [mmlon,idxlon] = min(closelon); %index of closest longitude grid to streamline point
                if latitude(jj) < 10 && latitude(jj) > -10 %use first-mode equatorial rossby phase speed for near equator
                    rossbyspeed(jj) = phspeedfine(idxlat,idxlon)/3;
                else %outside equatorial region
                    %calculate beta and f at this latitude
                    beta = 2*7.292*10^(-5)*cosd(latitude(jj))/(6371*10^3);
                    f = 2*7.292*10^(-5)*sind(latitude(jj));
                    rossbyspeed(jj) = beta*(phspeedfine(idxlat,idxlon)/f)^2;
                end
                
                if jj == length(ostrme(ii).lat)
                    rossbyspeed(jj) = rossbyspeed(jj-1);
                    rossbyspeed_plusU(jj) = rossbyspeed(jj-1) - ostrme(ii).u(jj-1);
                else
                    %add background u velocity
                    rossbyspeed_plusU(jj) = rossbyspeed(jj) - ostrme(ii).u(jj);
                    %project onto west direction
                    [~,phaseangle] = sw_dist([latitude(jj),latitude(jj+1)],[longitude(jj),longitude(jj+1)]);
                    rossbyspeed(jj) = rossbyspeed(jj)*cosd(phaseangle+180); 
                    rossbyspeed_plusU(jj) = -rossbyspeed_plusU(jj)*cosd(phaseangle);
                end
            end
            ostrme(ii).rossbyspeed = rossbyspeed;
            ostrme(ii).rossbyspeed_plusU = rossbyspeed_plusU;
        end
        save(strmnames{kk},'ostrme');
    end
    clear ostrm
    clear ostrme
end