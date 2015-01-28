clear all
close all
clc

path(path,'/home/sbrus/Codes/SomeDGMaterials')
path(path,'/home/sbrus/Codes/SomeDGMaterials/devscript')
path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

% direc = '/home/sbrus/Codes/dgswe/grids/';
% name = 'inlet1.grd';

% direc = '/home/sbrus/Codes/rimls/';
% name = 'galveston_tri.grd';

direc = '/home/sbrus/Codes/rimls/';
name = 'beaufort.grd';


[ect,xy,hb,~,~,~,~] = readfort14([direc,name]);

nn = length(xy);

fid = fopen('edge_nodes.d');

N = fscanf(fid,'%g %g',2) ;
nends = N(1)*N(2);
ends = fscanf(fid,'%g %g %g \n', [3 nends])' ;

fid4 = fopen('interior_nodes.d');

N = fscanf(fid4,'%g %g',2) ;
ninds = N(1)*N(2);
inds = fscanf(fid4,'%g %g %g \n', [3 ninds])' ;

fid2 = fopen('centers.d');

ne = fscanf(fid2,'%g',1) ;
cnds = fscanf(fid2,'%g %g %g \n', [3 ne])' ;

fid3 = fopen('normals.d');

ne = fscanf(fid3,'%g',1) ;
nrm = fscanf(fid3,'%g %g %g \n', [3 ne])' ;

fid5 = fopen('rimls_edge_nodes.d');

N = fscanf(fid5,'%g %g',2) ;
nends = N(1)*N(2);
rimls_ends = fscanf(fid5,'%g %g %g \n', [3 nends])' ;

fid6 = fopen('rimls_interior_nodes.d');

N = fscanf(fid6,'%g %g',2) ;
ninds = N(1)*N(2);
rimls_inds = fscanf(fid6,'%g %g %g \n', [3 ninds])' ;

fid7 = fopen('rimls_vertex_nodes.d');

N = fscanf(fid7,'%g %g',2) ;
ninds = N(1)*N(2);
rimls_vnds = fscanf(fid7,'%g %g %g \n', [3 ninds])' ;

xpts = vertcat(xy(:,1),ends(:,1),inds(:,1));
ypts = vertcat(xy(:,2),ends(:,2),inds(:,2));
hpts = vertcat(hb,ends(:,3),inds(:,3));
% xpts = vertcat(ends(:,1),inds(:,1));
% ypts = vertcat(ends(:,2),inds(:,2));
% hpts = vertcat(ends(:,3),inds(:,3));

rimls_xpts = vertcat(rimls_vnds(:,1),rimls_ends(:,1),rimls_inds(:,1));
rimls_ypts = vertcat(rimls_vnds(:,2),rimls_ends(:,2),rimls_inds(:,2));
rimls_hpts = vertcat(rimls_vnds(:,3),rimls_ends(:,3),rimls_inds(:,3));

elem = 'off';
node = 'off';
lcolor = lines;

DEToV = zeros(ne,4) ;
DEToV(:,1) = 3 ;
DEToV(:,2:4) = ect;

% figure(1);
% hold on
% drawNGonMesh4( xy, DEToV, 'ElNum', elem, 'NodeNum', node, lcolor(1,:) );

% plot(xy(:,1),xy(:,2),'o')
% plot(ends(:,1),ends(:,2),'ro')
% plot(inds(:,1),inds(:,2),'ro')
% plot(cnds(:,1),cnds(:,2),'b.')

% figure(5)
% pdeplot( xy', [], ect', 'xydata',hb,'zdata',hb, 'colormap', 'jet', 'mesh','on') ;

figure(2)
ect2 = delaunay(xpts,ypts);
nodes = [xpts ypts];
pdeplot( nodes', [], ect2', 'xydata',hpts,'zdata',hpts, 'colormap', 'jet', 'mesh','off') ;
% zlim([0 20])
% caxis([0 20])
% xlim([-77.5 -75])
% ylim([34.5 36.5])
% set(gca,'CameraPosition',[-82.18446985898792 33.028568423134736 174.72782070926633])
% set(gca,'CameraTarget',[-76.25 35.5 10.0])

% figure(3)
% quiver3(cnds(:,1),cnds(:,2),cnds(:,3),nrm(:,1),nrm(:,2),nrm(:,3))


figure(4)
ect3 = delaunay(rimls_xpts,rimls_ypts);
nodes = [rimls_xpts rimls_ypts];
pdeplot( nodes', [], ect3', 'xydata',rimls_hpts,'zdata',rimls_hpts, 'colormap', 'jet', 'mesh','off') ;
% zlim([0 20])
% caxis([0 20])
% xlim([-77.5 -75])
% ylim([34.5 36.5])
% set(gca,'CameraPosition',[-82.18446985898792 33.028568423134736 174.72782070926633])
% set(gca,'CameraTarget',[-76.25 35.5 10.0])





rmpath '/home/sbrus/Codes/SomeDGMaterials'
rmpath '/home/sbrus/Codes/SomeDGMaterials/devscript'
rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'
