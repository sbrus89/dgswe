close all
clear all
clc

% nodes_file = 'coarse.grd';
nodes_file = 'coarse_x2.grd';

fid = fopen(nodes_file) ;
agrid = fgetl(fid) ;
N = fscanf(fid,'%g %g',2) ;
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;
iv = sort(Val(:,1)) ;
Val = Val(iv,:) ; 
xy = Val(:,2:3);

idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;
iv = sort(idx(:,1)) ;
idx = idx(iv,:) ; 
ect = idx(:,3:end);




% grid_file = '/home/sbrus/data-drive/galveston_SL18/grid_dev/v17_cart/galveston_SL18_cart.grd';
grid_file = 'coarse_relaxed.grd';

fid = fopen(grid_file) ;
agrid = fgetl(fid) ;
N = fscanf(fid,'%g %g',2) ;
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;
iv = sort(Val(:,1)) ;
Val = Val(iv,:) ; 
xy_grid = Val(:,2:3);

idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;
iv = sort(idx(:,1)) ;
idx = idx(iv,:) ; 
ect_grid = idx(:,3:end);

hold all
% pdeplot( xy_grid', [], ect_grid','mesh','on') ;
plot(xy_grid(:,1),xy_grid(:,2),'b.')
pdeplot( xy', [], ect','mesh','on') ;
plot(xy(:,1),xy(:,2),'r.')
axis equal