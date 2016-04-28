clc
close all
clear all

% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri/p1/';};
%          
% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri/p2/';};
         
% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p3/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p3/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p3/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri/p3/';};         
%  labels = {'adcirc','adcircx64','quad spline channel','quad spline','quad','tri'};

% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p3/';};
         
% sta_direc = {%'/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p1/';
%          %    '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p3/';};
%          
% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p3/';};         
          
% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/tri/p1/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri/p2/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri/p3/';};

sta_direc = {%'/home/sbrus/data-drive/galveston/dgswe/tri_spline/p1/';
            % '/home/sbrus/data-drive/galveston/dgswe/tri_spline/p2/';
             '/home/sbrus/data-drive/galveston/dgswe/tri_spline/p3/';};

labels = {'adcirc';
          'adcircx64';
          %'DG p1';
          %'DG p2';
          'DG p3'}; 

% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/tri/p3/';
%              '/home/sbrus/data-drive/galveston/dgswe/tri_spline/p3/';};
       
% sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2/M2/p3/';
%              '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/p3/';};         
%          
% labels = {'adcirc','adcircx64','p3 straight','p3 spline'}; 
         
adc_direc1 = '/home/sbrus/data-drive/galveston/adcirc/original/ESL0/';
adc_direc2 = '/home/sbrus/data-drive/galveston/adcirc/refinedx64/ESL0/';
        
% grid_name = '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p1/galveston2_plot.grd';
grid_name = '/home/sbrus/data-drive/galveston/dgswe/tri/p1/galveston_tri.grd';

% sta = 1:401
sta = 1:81; %middle 
% sta = 82:163; %left
%  sta = 164:246; %right

nsnap = 24;
nruns = length(sta_direc);

figure(1)
hold all
figure(2)
hold all

plot_adcirc_stations(adc_direc1,nsnap,sta)
plot_adcirc_stations(adc_direc2,nsnap,sta)

for run = 1:nruns
fid_H = fopen([sta_direc{run},'station_H.d']);
fid_Qx = fopen([sta_direc{run},'station_Qx.d']);
fid_Qy = fopen([sta_direc{run},'station_Qy.d']);

nsta = fscanf(fid_H,' %g ', 1);
nsta = fscanf(fid_Qx,' %g ', 1);
nsta = fscanf(fid_Qy,' %g ', 1);

H = zeros(nsta,nsnap);
hb = zeros(nsta,nsnap);
Qx = zeros(nsta,nsnap);
Qy = zeros(nsta,nsnap);

xyH = zeros(nsta,3);
xyQx = zeros(nsta,3);
xyQy = zeros(nsta,3);

snap = 0;
while ~feof(fid_H) && snap < nsnap
    
    snap = snap + 1;
    
    th = fscanf(fid_H,' %g ', 1);
    xyH = fscanf(fid_H,' %g ', [4 nsta])'; 
    
    tqx = fscanf(fid_Qx,' %g ', 1); 
    xyQx = fscanf(fid_Qx,' %g ', [3 nsta])'; 
    
    t(snap) = fscanf(fid_Qy,' %g ', 1); 
    xyQy = fscanf(fid_Qy,' %g ', [3 nsta])';
    
    H(:,snap) = xyH(:,3);
    hb(:,snap) = xyH(:,4);
    Qx(:,snap) = xyQx(:,3);
    Qy(:,snap) = xyQy(:,3);
end

u = Qx./H;
v = Qy./H;

% nsnap = snap;

for snap = nsnap
    figure(1)
    plot(H(sta,snap)-hb(sta,snap));
    figure(2)
    plot(sqrt(u(sta,snap).^2 + v(sta,snap).^2))

end

end

figure(1)
% legend(labels,'Location','Southwest','FontSize',10)
xlabel('station')
ylabel('surface elevation (m)')
set(gcf, 'units', 'inches', 'pos', [0 0 3.5 2.8])
set(gcf,'PaperPositionMode','auto') ;
print('-dpdf',['~/Presentations/proposal/figures/galveston_zeta_',num2str(nsnap),'.pdf'])
figure(2)
% legend(labels,'Location','Northwest','FontSize',10)
xlabel('station')
ylabel('velocity (m/s)')
set(gcf, 'units', 'inches', 'pos', [0 0 3.5 2.8])
set(gcf,'PaperPositionMode','auto') ;
print('-dpdf',['~/Presentations/proposal/figures/galveston_vel_',num2str(nsnap),'.pdf'])

% figure
% hold on
% plot_grid(grid_name);
% plot(xyH(sta,1),xyH(sta,2),'ro','MarkerFaceColor','r','MarkerSize',2)
% % quiver(xyH(sta,1),xyH(sta,2),u(sta,snap),v(sta,snap),'b')
