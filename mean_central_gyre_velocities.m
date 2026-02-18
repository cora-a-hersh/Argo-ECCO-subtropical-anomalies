%% load in climatologies
addpath /atalanta/home/swijffels/work/argo/gridNSF/climfiles
addpath /atalanta/home/chersh/toolbox
addpath /atalanta/home/chersh/toolbox/cmocean
addpath /atalanta/home/swijffels/work/argo/
addpath /atalanta/home/swijffels/toolbox/seawater
addpath /atalanta/home/swijffels/toolbox/csirolib
addpath /atalanta/home/swijffels/work/argo/matlab
addpath /atalanta/home/swijffels/toolbox/subaxis
addpath /atalanta/home/swijffels/toolbox/susan

%a = load('/home/swijffels/work/argo/gridNSF/climfiles/clim_pres_geovel.mat');
a = load('/atalanta/home/swijffels/work/argo/gridNSF/geovel/climatology_mean_geovel_nondivref_cora_argo_sla_ug_lscovtar_huber25_200.mat');
e = load('/atalanta/home/swijffels/work/argo/gridNSF/ecco4/climatology_mean_geovel_1000dbref_iter129_newgrid.mat');
e_on_agrid = load('/atalanta/home/chersh/SpiceAnomalies/compare_argo_ecco_mean_properties/eccoclim_on_argo_grid.mat');

%% get lat/lon/sig grids
e_xi = e.xi;
e_yi = e.yi;
a_xi = a.xi;
a_yi = a.yi;
sig1grid = a.sig1grid;
[a_XI,a_YI] = meshgrid(a_xi,a_yi);
[e_XI,e_YI] = meshgrid(e_xi,e_yi);
pgrid = a.pgrid;

%% mean velocity/depth profile in center of basin

%North Pacific
lonrange = [190,210];
latrange = [15,30];
%lonrange = [330,340];
%latrange = [-30,-15];

smoothing = false;

[~,ilonstart] = min(abs((a_xi - lonrange(1))));
[~,ilonend] = min(abs((a_xi - lonrange(2))));
[~,ilatstart] = min(abs((a_yi - latrange(1))));
[~,ilatend] = min(abs((a_yi - latrange(2))));

uga = squeeze(a.ug(ilonstart:ilonend,ilatstart:ilatend,:));
vga = squeeze(a.vg(ilonstart:ilonend,ilatstart:ilatend,:));

if smoothing == true
    uga = imgaussfilt(uga,5);
    vga = imgaussfilt(vga,5);
end

ugamean = squeeze(mean(uga,1,'omitnan'));
ugamean = squeeze(mean(ugamean,1,'omitnan'));
vgamean = squeeze(mean(vga,1,'omitnan'));
vgamean = squeeze(mean(vgamean,1,'omitnan'));

%interpolated ecco velocities to argo x/y grid
ude = squeeze(e_on_agrid.ude_interp(ilonstart:ilonend,ilatstart:ilatend,:));
vde = squeeze(e_on_agrid.vde_interp(ilonstart:ilonend,ilatstart:ilatend,:));

udemean = squeeze(mean(ude,1,'omitnan'));
udemean = squeeze(mean(udemean,1,'omitnan'));
vdemean = squeeze(mean(vde,1,'omitnan'));
vdemean = squeeze(mean(vdemean,1,'omitnan'));

figure('Position',[10 10 900 800])
subaxis(2,2,1,'SH',0.04,'SV',0.03)
plot(ugamean,pgrid,'r','LineWidth',2)
hold on
plot(vgamean,pgrid,'b','LineWidth',2)
plot(udemean,pgrid,'r--','LineWidth',2)
plot(vdemean,pgrid,'b--','LineWidth',2)
plot([0 0],[pgrid(1) pgrid(end)],'k:')
set(gca,'XTick',[])
ylabel('dbar')
xlim([-0.03 0.03])
axis ij
vax = axis;
text(vax(1)+0.005,vax(4)-250,'North Pacific')

%North Atlantic
lonrange = [310,330];
latrange = [15,30];

smoothing = false;

[~,ilonstart] = min(abs((a_xi - lonrange(1))));
[~,ilonend] = min(abs((a_xi - lonrange(2))));
[~,ilatstart] = min(abs((a_yi - latrange(1))));
[~,ilatend] = min(abs((a_yi - latrange(2))));

uga = squeeze(a.ug(ilonstart:ilonend,ilatstart:ilatend,:));
vga = squeeze(a.vg(ilonstart:ilonend,ilatstart:ilatend,:));

if smoothing == true
    uga = imgaussfilt(uga,5);
    vga = imgaussfilt(vga,5);
end

ugamean = squeeze(mean(uga,1,'omitnan'));
ugamean = squeeze(mean(ugamean,1,'omitnan'));
vgamean = squeeze(mean(vga,1,'omitnan'));
vgamean = squeeze(mean(vgamean,1,'omitnan'));

%interpolated ecco velocities to argo x/y grid
ude = squeeze(e_on_agrid.ude_interp(ilonstart:ilonend,ilatstart:ilatend,:));
vde = squeeze(e_on_agrid.vde_interp(ilonstart:ilonend,ilatstart:ilatend,:));

udemean = squeeze(mean(ude,1,'omitnan'));
udemean = squeeze(mean(udemean,1,'omitnan'));
vdemean = squeeze(mean(vde,1,'omitnan'));
vdemean = squeeze(mean(vdemean,1,'omitnan'));

subaxis(2,2,2)
plot(ugamean,pgrid,'r','LineWidth',2)
hold on
plot(vgamean,pgrid,'b','LineWidth',2)
plot(udemean,pgrid,'r--','LineWidth',2)
plot(vdemean,pgrid,'b--','LineWidth',2)
plot([0 0],[pgrid(1) pgrid(end)],'k:')
set(gca,'XTick',[])
set(gca,'YTick',[])
xlim([-0.03 0.03])
axis ij
vax = axis;
text(vax(1)+0.005,vax(4)-250,'North Atlantic')

%South Pacific
lonrange = [200,220];
latrange = [-30,-15];

smoothing = false;

[~,ilonstart] = min(abs((a_xi - lonrange(1))));
[~,ilonend] = min(abs((a_xi - lonrange(2))));
[~,ilatstart] = min(abs((a_yi - latrange(1))));
[~,ilatend] = min(abs((a_yi - latrange(2))));

uga = squeeze(a.ug(ilonstart:ilonend,ilatstart:ilatend,:));
vga = squeeze(a.vg(ilonstart:ilonend,ilatstart:ilatend,:));

if smoothing == true
    uga = imgaussfilt(uga,5);
    vga = imgaussfilt(vga,5);
end

ugamean = squeeze(mean(uga,1,'omitnan'));
ugamean = squeeze(mean(ugamean,1,'omitnan'));
vgamean = squeeze(mean(vga,1,'omitnan'));
vgamean = squeeze(mean(vgamean,1,'omitnan'));

%interpolated ecco velocities to argo x/y grid
ude = squeeze(e_on_agrid.ude_interp(ilonstart:ilonend,ilatstart:ilatend,:));
vde = squeeze(e_on_agrid.vde_interp(ilonstart:ilonend,ilatstart:ilatend,:));

udemean = squeeze(mean(ude,1,'omitnan'));
udemean = squeeze(mean(udemean,1,'omitnan'));
vdemean = squeeze(mean(vde,1,'omitnan'));
vdemean = squeeze(mean(vdemean,1,'omitnan'));

subaxis(2,2,3)
plot(ugamean,pgrid,'r','LineWidth',2)
hold on
plot(vgamean,pgrid,'b','LineWidth',2)
plot(udemean,pgrid,'r--','LineWidth',2)
plot(vdemean,pgrid,'b--','LineWidth',2)
plot([0 0],[pgrid(1) pgrid(end)],'k:')
axis ij
xlim([-0.03 0.03])
xlabel('m/s')
ylabel('dbar')
vax = axis;
text(vax(1)+0.005,vax(4)-250,'South Pacific')

%South Indian
lonrange = [70,90];
latrange = [-30,-15];

smoothing = false;

[~,ilonstart] = min(abs((a_xi - lonrange(1))));
[~,ilonend] = min(abs((a_xi - lonrange(2))));
[~,ilatstart] = min(abs((a_yi - latrange(1))));
[~,ilatend] = min(abs((a_yi - latrange(2))));

uga = squeeze(a.ug(ilonstart:ilonend,ilatstart:ilatend,:));
vga = squeeze(a.vg(ilonstart:ilonend,ilatstart:ilatend,:));

if smoothing == true
    uga = imgaussfilt(uga,5);
    vga = imgaussfilt(vga,5);
end

ugamean = squeeze(mean(uga,1,'omitnan'));
ugamean = squeeze(mean(ugamean,1,'omitnan'));
vgamean = squeeze(mean(vga,1,'omitnan'));
vgamean = squeeze(mean(vgamean,1,'omitnan'));

%interpolated ecco velocities to argo x/y grid
ude = squeeze(e_on_agrid.ude_interp(ilonstart:ilonend,ilatstart:ilatend,:));
vde = squeeze(e_on_agrid.vde_interp(ilonstart:ilonend,ilatstart:ilatend,:));

udemean = squeeze(mean(ude,1,'omitnan'));
udemean = squeeze(mean(udemean,1,'omitnan'));
vdemean = squeeze(mean(vde,1,'omitnan'));
vdemean = squeeze(mean(vdemean,1,'omitnan'));

subaxis(2,2,4)
plot(ugamean,pgrid,'r','LineWidth',2)
hold on
plot(vgamean,pgrid,'b','LineWidth',2)
plot(udemean,pgrid,'r--','LineWidth',2)
plot(vdemean,pgrid,'b--','LineWidth',2)
plot([0 0],[pgrid(1) pgrid(end)],'k:')
set(gca,'YTick',[])
axis ij
xlim([-0.03 0.03])
xlabel('m/s')
legend('Argo u','Argo v','ECCO u','ECCO v')
vax = axis;
text(vax(1)+0.005,vax(4)-250,'South Indian')

%saveas(gcf,'/home/chersh/SpiceAnomalies/Paper1_figures/mean_central_gyre_velocities.png')
%saveas(gcf,'/home/chersh/SpiceAnomalies/Paper1_figures/mean_central_gyre_velocities.fig')
