close all
clear all
clc

% file =  '/home/sbrus/Codes/dgswe/grids/inlet1_ctp3.cb'
% file = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_quad_x4_ctp3.cb';
% file = '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p3/ctp3/hbp3/decomp48/PE0009/fort.cb';
file = '/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/spline_modified/coarse_x2_ctp3.cb';

fid = fopen(file);
header = 1;
while header == 1
    line = fgetl(fid);
    disp(line)
    if length(line) > 10 && strcmp(line(1:10),'!!!!!!!!!!');
        header = 0;
    end
end

nbou = fscanf(fid,'%d', 1); fgetl(fid);
n = fscanf(fid,'%d %d', 2); fgetl(fid);
nvel = n(1);
ctp = n(2);

figure(5)
hold on

for bou = 1:nbou
    
    n = fscanf(fid,'%d %d', 2); fgetl(fid);
    disp(n)
    nnds = n(1);
    btype = n(2);
    
    if nnds <= 0
        continue
    end
    
    if ctp == 2
        nodes = fscanf(fid,'%d %g %g ', [ctp+1,2*nnds-2])';
    elseif ctp == 3
        nodes = fscanf(fid,'%d %g %g %g', [ctp+1,2*nnds-2])';    
    end
    disp(nodes)
    last = fscanf(fid,'%d %g', [2,2])';
    
    xtmp = nodes(1:2:2*nnds-2,3:end);
    xnodes = vertcat(reshape(xtmp,[(ctp-1)*(nnds-1),1]));
    xtmp = nodes(1:2:2*nnds-2,2);
    xvert =  vertcat(reshape(xtmp,[(nnds-1),1]),last(1,2));
    
    ytmp = nodes(2:2:2*nnds-1,3:end);
    ynodes = vertcat(reshape(ytmp,[(ctp-1)*(nnds-1),1]));  
    ytmp = nodes(2:2:2*nnds-1,2);
    yvert =  vertcat(reshape(ytmp,[(nnds-1),1]),last(2,2));    
    
    plot(xnodes,ynodes,'kx')
    plot(xvert,yvert,'ko')
end

axis equal

fclose(fid);
