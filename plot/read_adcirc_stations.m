function [z,vel,time] = read_adcirc_stations(direc)



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
snap = 0;
while ~feof(f62) && snap < nsnap
   snap = snap + 1;
   
   t = fscanf(f62,'%f %d \n',2);
   time(snap) = t(1);
   
   val = fscanf(f62, '%d %f %f\n', [3, nsta])';
   u(:,snap) = val(:,2);
   v(:,snap) = val(:,3);
end

vel = sqrt(u.^2 + v.^2);
fclose(f62);
