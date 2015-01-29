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
% name = 'beaufort.grd';
% name = 'galveston_tri.grd';
% name = 'shin.grd';

direc = '/home/sbrus/data-drive/EC2001/';
name = 'fort.14';

mesh = 'on';

% % entire area
% xbox = [min(xpts) max(xpts)];
% ybox = [min(ypts) max(ypts)];

% % beaufort detail
% xbox = [-77 -75.3];
% ybox = [35.7 36.5];

% %EC2001 Bermuda
% xbox = [-66.08 -63.38];
% ybox = [30.76 34.04];

% %EC2001 PR
% xbox = [-68.09 -64.79];
% ybox = [16.92 19.04];

%EC2001 Cayman
xbox = [-86.52 -79.97];
ybox = [16.75 22.99];

% xbox = [-74.87 -73.09];
% ybox = [35.14 37.70];


[ect,xy,hb,nelnds,opedat,boudat,title,bvnds] = readfort14([direc,name]);

nn = length(xy);

fid = fopen('edge_nodes.d');
N = fscanf(fid,'%g %g',2) ;
nends = N(1)*N(2);
ends = fscanf(fid,'%g %g %g \n', [3 nends])' ;
fclose(fid);

fid4 = fopen('interior_nodes.d');
N = fscanf(fid4,'%g %g',2) ;
ninds = N(1)*N(2);
inds = fscanf(fid4,'%g %g %g \n', [3 ninds])' ;
fclose(fid4);

fid2 = fopen('centers.d');
ne = fscanf(fid2,'%g',1) ;
cnds = fscanf(fid2,'%g %g %g \n', [3 ne])' ;
fclose(fid2);

fid3 = fopen('normals.d');
ne = fscanf(fid3,'%g',1) ;
nrm = fscanf(fid3,'%g %g %g \n', [3 ne])' ;
fclose(fid3);

fid8 = fopen('boundary_nodes.d');
N = fscanf(fid8,'%g %g',2) ;
nends = N(1)*N(2);
bends = fscanf(fid8,'%g \n', [1 nends])' ;
fclose(fid8);


fid5 = fopen('rimls_edge_nodes.d');
N = fscanf(fid5,'%g %g',2) ;
nends = N(1)*N(2);
rimls_ends = fscanf(fid5,'%g %g %g \n', [3 nends])' ;
fclose(fid5);

fid6 = fopen('rimls_interior_nodes.d');
N = fscanf(fid6,'%g %g',2) ;
ninds = N(1)*N(2);
rimls_inds = fscanf(fid6,'%g %g %g \n', [3 ninds])' ;
fclose(fid6);

fid7 = fopen('rimls_vertex_nodes.d');
N = fscanf(fid7,'%g %g',2) ;
ninds = N(1)*N(2);
rimls_vnds = fscanf(fid7,'%g %g %g \n', [3 ninds])' ;
fclose(fid7);




xpts = vertcat(xy(:,1),ends(:,1),inds(:,1));
ypts = vertcat(xy(:,2),ends(:,2),inds(:,2));
hpts = vertcat(hb,ends(:,3),inds(:,3));

rimls_xpts = vertcat(rimls_vnds(:,1),rimls_ends(:,1),rimls_inds(:,1));
rimls_ypts = vertcat(rimls_vnds(:,2),rimls_ends(:,2),rimls_inds(:,2));
rimls_hpts = vertcat(rimls_vnds(:,3),rimls_ends(:,3),rimls_inds(:,3));

bnds_flag = vertcat(bvnds,bends,0*inds(:,1));




in = inpolygon (xpts,ypts,xbox,ybox);

xplot = xpts(in);
yplot = ypts(in);
hplot = hpts(in);
bnds_plot = bnds_flag(in);

in = inpolygon (rimls_xpts,rimls_ypts,xbox,ybox);

rimls_xplot = rimls_xpts(in);
rimls_yplot = rimls_ypts(in);
rimls_hplot = rimls_hpts(in);
rimls_bnds_plot = bnds_flag(in);

in = inpolygon (xy(:,1),xy(:,2),xbox,ybox);

xvplot = xy(in,1);
yvplot = xy(in,2);
hvplot = hb(in);
bvnds_plot = bvnds(in);




figure(1)
pdeplot( xy', [], ect', 'xydata',hb, 'colormap', 'jet', 'mesh','on') ;
hold on 
plot([xbox(1) xbox(2)],[ybox(1) ybox(1)],'r','LineWidth',5)
plot([xbox(1) xbox(2)],[ybox(2) ybox(2) ],'r','LineWidth',5)
plot([xbox(1) xbox(1)],[ybox(1) ybox(2) ],'r','LineWidth',5)
plot([xbox(2) xbox(2)],[ybox(1) ybox(2) ],'r','LineWidth',5)
axis image
grid on

figure(2)
ect1 = delaunay(xvplot,yvplot);
ect1c = clean_ect(ect1,bvnds_plot);
nodes = [xvplot yvplot];
pdeplot( nodes', [], ect1c', 'xydata',hvplot, 'zdata',hvplot,'colormap', 'jet', 'mesh','on') ;
grid on

figure(3)
ect2 = delaunay(xplot,yplot);
ect2c = clean_ect(ect2,bnds_plot);
nodes = [xplot yplot];
pdeplot( nodes', [], ect2c', 'xydata',hplot,'zdata',hplot, 'colormap', 'jet', 'mesh',mesh) ;
grid on


figure(4)
ect3 = delaunay(rimls_xplot,rimls_yplot);
ect3c = clean_ect(ect3,rimls_bnds_plot);
nodes = [rimls_xplot rimls_yplot];
pdeplot( nodes', [], ect3c', 'xydata',rimls_hplot,'zdata',rimls_hplot, 'colormap', 'jet', 'mesh',mesh) ;
grid on

% figure(5)
% quiver3(cnds(:,1),cnds(:,2),cnds(:,3),nrm(:,1),nrm(:,2),nrm(:,3))

figure(6);
pdeplot( [xvplot yvplot]', [], ect1c','mesh','on') ;
axis image

hold on
% plot(xy(:,1),xy(:,2),'o')
% plot(ends(:,1),ends(:,2),'ro')
% plot(inds(:,1),inds(:,2),'ro')
% plot(cnds(:,1),cnds(:,2),'b.')
plot(xplot,yplot,'ro')


disp('Max absolute difference between rimls and  linear interpolation')
disp(max(abs(rimls_hpts-hpts)))
disp('Max relative difference between rimls and linear interpolation')
disp(max(abs(rimls_hpts-hpts))/max(hpts)*100)



rmpath '/home/sbrus/Codes/SomeDGMaterials'
rmpath '/home/sbrus/Codes/SomeDGMaterials/devscript'
rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'
