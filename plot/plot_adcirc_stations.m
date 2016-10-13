function plot_adcirc_stations(direc,tsnap,sta,ts_sta,color)

% close all
% clear all
% clc

% direc = '/home/sbrus/data-drive/galveston/adcirc/refinedx64/ESL0/';
% direc = '/home/sbrus/data-drive/galveston/adcirc/original/ESL0/';

f61 = fopen([direc,'fort.61']);

line = fgetl(f61);

n = fscanf(f61,'%d %d %f %d %d %s %d \n',7);

nsnap = n(1,1);
nsta = n(2,1);

z = zeros(nsta,nsnap);
snap = 0;
while ~feof(f61) && snap < nsnap
   snap = snap + 1;
   t = fscanf(f61,'%f %d \n',2);
   time = t(1);
   
   val = fscanf(f61, '%d %f \n', [2, nsta])';
   z(:,snap) = val(:,2);
end

fclose(f61);



f62 = fopen([direc,'fort.62']);
line = fgetl(f62);

n = fscanf(f62,'%d %d %f %d %d %s %d \n',7);

nsnap = n(1,1);
nsta = n(2,1);

u = zeros(nsta,nsnap);
v = zeros(nsta,nsnap);
for snap = 1:nsnap
   t = fscanf(f62,'%f %d \n',2);
   time(snap) = t(1);
   
   val = fscanf(f62, '%d %f %f\n', [3, nsta])';
   u(:,snap) = val(:,2);
   v(:,snap) = val(:,3);
end

vel = sqrt(u.^2 + v.^2);
fclose(f62);


for snap = tsnap-1
    figure(1)
    plot(z(sta,snap),'Color',color)
    figure(2)
    plot(vel(sta,snap),'Color',color)
end


figure(3)
plot(time,z(ts_sta,:),'Color',color)
figure(4)
plot(time,vel(ts_sta,:),'Color',color)