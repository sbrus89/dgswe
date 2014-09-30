clc
% close all
clear all

% sta_direc = {'/home/sbrus/Codes/stations/'};
sta_direc = {'/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/';
             '/home/sbrus/data-drive/galveston/dgswe/quad2_spline/';
             '/home/sbrus/data-drive/galveston/dgswe/quad2/M2/';
             '/home/sbrus/data-drive/galveston/dgswe/tri/';};
grid_name = '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/galveston2_plot.grd';

nruns = length(sta_direc);

for run = 1:nruns
fid_H = fopen([sta_direc{run},'station_H.d']);
fid_Qx = fopen([sta_direc{run},'station_Qx.d']);
fid_Qy = fopen([sta_direc{run},'station_Qy.d']);

nsta = fscanf(fid_H,' %g ', 1);
nsta = fscanf(fid_Qx,' %g ', 1);
nsta = fscanf(fid_Qy,' %g ', 1);

nsnap = 50;

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

nsnap = snap;
sta = 1:81;
% sta = 82:163;
% sta = 164:246;
for snap = 46
    figure(1)
    plot(H(sta,snap)-hb(sta,snap));
    figure(2)
    plot(sqrt(u(sta,snap).^2 + v(sta,snap).^2))
end
figure(1)
hold all
figure(2)
hold all
end

figure
hold on
plot_grid(grid_name);
plot(xyH(sta,1),xyH(sta,2),'ro')
quiver(xyH(sta,1),xyH(sta,2),u(sta,snap),v(sta,snap),'b')
