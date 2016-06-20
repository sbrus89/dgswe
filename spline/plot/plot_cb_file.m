close all
clear all
clc

file = '/home/sbrus/data-drive/galveston_spline/grids/galveston_tri_x64_ctp3.cb';
fid = fopen(file);
header = fgetl(fid);

nbou = fscanf(fid,'%d', 1); fgetl(fid);
n = fscanf(fid,'%d %d', 2); fgetl(fid);
nvel = n(1);
ctp = n(2);

figure
hold on

for bou = 1:nbou
    
    n = fscanf(fid,'%d %d', 2); fgetl(fid);
    nnds = n(1);
    btype = n(2);
    
    nodes = fscanf(fid,'%d %g %g %g', [ctp+1,2*nnds-2])';
    last = fscanf(fid,'%d %g', [2,2])';
    
    xtmp = nodes(1:2:2*nnds-2,3:end);
    xnodes = vertcat(reshape(xtmp,[(ctp-1)*(nnds-1),1]));
    xtmp = nodes(1:2:2*nnds-2,2);
    xvert =  vertcat(reshape(xtmp,[(nnds-1),1]),last(1,2));
    
    ytmp = nodes(2:2:2*nnds-1,3:end);
    ynodes = vertcat(reshape(ytmp,[(ctp-1)*(nnds-1),1]));  
    ytmp = nodes(2:2:2*nnds-1,2);
    yvert =  vertcat(reshape(ytmp,[(nnds-1),1]),last(2,2));    
    
    plot(xnodes,ynodes,'x')
    plot(xvert,yvert,'o')
end

axis equal

fclose(fid);
