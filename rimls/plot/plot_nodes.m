clear all
close all
clc

% path(path,'/home/sbrus/Codes/SomeDGMaterials')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/devscript')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

out_dir = '../output/';

finp = fopen('../work/rimls.inp');

n = 0;
while n < 5 && ~feof(finp)
   temp = strtrim(fgetl(finp));
     
   if isempty(temp) || temp(1) == '!'
       
   else
      n = n + 1;
      temp2 = strtrim(strsplit(temp));
      names{n} = temp2{1};
   end
end

base_name = names{1};
eval_name = names{5};


mesh = 'on';

% entire area
xbox = [0 0];
ybox = [0 0];

% galveston trinity bay
% xbox = [3.152e5 3.4e5];
% ybox = [3.272e6 3.3e6];


% converging channel constriction
% xbox = [1875 2125];
% ybox = [50 160];

% % galveston SL18+TX33 inlet
% xbox = [-94.835 -94.667];
% ybox = [29.298 29.404];

% beaufort deep
% xbox = [-75.5 -73];
% ybox = [33 37];

% beaufort detail
% xbox = [-77 -75.3];
% ybox = [35.7 36.5];

% % beaufort channel
% xbox = [-76.8 -76.6];
% ybox = [34.7 35];

% %EC2001 Bermuda
% xbox = [-66.08 -63.38];
% ybox = [30.76 34.04];

% %EC2001 PR
% xbox = [-68.09 -64.79];
% ybox = [16.92 19.04];

% %EC2001 Cayman
% xbox = [-86.52 -79.97];
% ybox = [16.75 22.99];

% xbox = [-74.87 -73.09];
% ybox = [35.14 37.70];


[base_ect,base_xy,base_hb,base_nelnds,opedat,boudat,header,base_bvnds] = readfort14(base_name);
if strcmp(base_name,eval_name)
    eval_ect = base_ect;
    eval_xy = base_xy;
    eval_nelnds = base_nelnds;
    eval_hb = base_hb;
    eval_bvnds = base_bvnds;
else
    [eval_ect,eval_xy,eval_hb,eval_nelnds,opedat,boudat,header,eval_bvnds] = readfort14(eval_name);
end


fid = fopen([out_dir,'edge_nodes.d']);
nends = fscanf(fid,'%g',1) ;
ends = fscanf(fid,'%g %g %g \n', [3 nends])' ;
fclose(fid);

fid4 = fopen([out_dir,'interior_nodes.d']);
ninds = fscanf(fid4,'%g',1) ;
inds = fscanf(fid4,'%g %g %g \n', [3 ninds])' ;
fclose(fid4);

fid2 = fopen([out_dir,'base_nodes.d']);
nbnds = fscanf(fid2,'%g',1) ;
base_nds = fscanf(fid2,'%g %g %g \n', [3 nbnds])' ;
fclose(fid2);

fid3 = fopen([out_dir,'normals.d']);
ne = fscanf(fid3,'%g',1) ;
nrm = fscanf(fid3,'%g %g %g \n', [3 ne])' ;
fclose(fid3);

fid3 = fopen([out_dir,'base_centers.d']);
ne = fscanf(fid3,'%g',1) ;
base_centers = fscanf(fid3,'%g %g \n', [2 ne])' ;
fclose(fid3);

fid8 = fopen([out_dir,'boundary_nodes.d']);
nends = fscanf(fid8,'%g',1) ;
bends = fscanf(fid8,'%g \n', [1 nends])' ;
fclose(fid8);


fid5 = fopen([out_dir,'rimls_edge_nodes.d']);
nends = fscanf(fid5,'%g',1) ;
rimls_ends = fscanf(fid5,'%g %g %g \n', [3 nends])' ;
fclose(fid5);

fid6 = fopen([out_dir,'rimls_interior_nodes.d']);
ninds = fscanf(fid6,'%g',1) ;
rimls_inds = fscanf(fid6,'%g %g %g \n', [3 ninds])' ;
fclose(fid6);

fid7 = fopen([out_dir,'rimls_vertex_nodes.d']);
ninds = fscanf(fid7,'%g',1) ;
rimls_vnds = fscanf(fid7,'%g %g %g \n', [3 ninds])' ;
fclose(fid7);

fid9 = fopen([out_dir,'curved_element_nodes.d']);
nends = fscanf(fid9,'%g',1) ;
cel_nds = fscanf(fid9,'%g %g \n', [2 nends])' ;
fclose(fid9);


% fid8 = fopen([out_dir,'random_pts.d']);
% ne = fscanf(fid8,'%g',1) ;
% random_pts = fscanf(fid8,'%g %g %g \n', [3 ne])' ;
% fclose(fid8);




xpts = vertcat(eval_xy(:,1),ends(:,1),inds(:,1));
ypts = vertcat(eval_xy(:,2),ends(:,2),inds(:,2));
hpts = vertcat(eval_hb,ends(:,3),inds(:,3));

rimls_xpts = vertcat(rimls_vnds(:,1),rimls_ends(:,1),rimls_inds(:,1));
rimls_ypts = vertcat(rimls_vnds(:,2),rimls_ends(:,2),rimls_inds(:,2));
rimls_hpts = vertcat(rimls_vnds(:,3),rimls_ends(:,3),rimls_inds(:,3));

bnds_flag = vertcat(eval_bvnds,bends,0*inds(:,1));


if all(xbox) == 0 && all(ybox) == 0
    xbox = [min(xpts) max(xpts)];
    ybox = [min(ypts) max(ypts)];
end


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

in = inpolygon (base_xy(:,1),base_xy(:,2),xbox,ybox);

base_xvert_plot = base_xy(in,1);
base_yvert_plot = base_xy(in,2);
base_hvert_plot = base_hb(in);
base_bouvert_plot = base_bvnds(in);
[base_ect_plot,base_nelnds_plot] = filter_ect(base_ect,base_nelnds,in);

in = inpolygon (eval_xy(:,1),eval_xy(:,2),xbox,ybox);

eval_xvert_plot = eval_xy(in,1);
eval_yvert_plot = eval_xy(in,2);
[eval_ect_plot,eval_nelnds_plot] = filter_ect(eval_ect,eval_nelnds,in);

in = inpolygon(base_nds(:,1),base_nds(:,2),xbox,ybox);
base_xnds_plot = base_nds(in,1);
base_ynds_plot = base_nds(in,2);

in = inpolygon(cel_nds(:,1),cel_nds(:,2),xbox,ybox);
cel_xnds_plot = cel_nds(in,1);
cel_ynds_plot = cel_nds(in,2);

in = inpolygon(base_centers(:,1),base_centers(:,2),xbox,ybox);
base_xcenters_plot = base_centers(in,1);
base_ycenters_plot = base_centers(in,2);



% figure(99)
% plot(random_pts(:,1),random_pts(:,2),'b.')
% 
% figure(100);
% ect_rand = delaunay(random_pts(:,1),random_pts(:,2));
% nodes = [random_pts(:,1) random_pts(:,2)];
% pdeplot( nodes', [], ect_rand', 'xydata',-random_pts(:,3), 'zdata',-random_pts(:,3),'colormap', 'jet', 'mesh','on') ;




% Plot of full grid batymetry and zoom box
% figure(1)
% pdeplot( base_xy', [], base_ect', 'xydata',base_hb, 'colormap', 'jet', 'mesh','on') ;
% hold on 
% plot([xbox(1) xbox(2)],[ybox(1) ybox(1)],'r','LineWidth',5)
% plot([xbox(1) xbox(2)],[ybox(2) ybox(2) ],'r','LineWidth',5)
% plot([xbox(1) xbox(1)],[ybox(1) ybox(2) ],'r','LineWidth',5)
% plot([xbox(2) xbox(2)],[ybox(1) ybox(2) ],'r','LineWidth',5)
% xlabel('longitude')
% ylabel('latitude')
% axis image
% grid on


%%
% % Plot of grid and extra nodes in zoom box
% figure(2);
% hold on
% drawNGonMesh4(horzcat(base_xvert_plot,base_yvert_plot), horzcat(base_nelnds_plot,base_ect_plot), 'ElNum', 'off', 'NodeNum', 'off', 'b',1 )
% plot(base_xnds_plot,base_ynds_plot,'g.')
% % plot(base_xcenters_plot,base_ycenters_plot,'bo','MarkerSize',4,'MarkerFaceColor','b')
% drawNGonMesh4(horzcat(eval_xvert_plot,eval_yvert_plot), horzcat(eval_nelnds_plot,eval_ect_plot), 'ElNum', 'off', 'NodeNum', 'off', 'k',3 )
% plot(xplot,yplot,'ro','MarkerSize',9,'MarkerFaceColor','r')
% plot(cel_xnds_plot,cel_ynds_plot,'ko','MarkerSize',7,'MarkerFaceColor','k')
% xlabel('longitude')
% ylabel('latitude')
% title('grid and evaluation points')
% axis equal

% break

% Plot of base grid bathymetry nodes
figure(3)
ect2 = delaunay(base_nds(:,1),base_nds(:,2));
nodes = [base_nds(:,1),base_nds(:,2)];
pdeplot( nodes', [], ect2', 'xydata',base_nds(:,3), 'zdata',base_nds(:,3),'colormap', 'jet', 'mesh','on') ;
xlabel('longitude')
ylabel('latitude')
zlabel('bathymetry (m)')
title('base grid bathymetry')
grid on

%%
% Plot of linearly interpolated nodes
figure(4)
ect2 = delaunay(xplot,yplot);
ect2c = clean_ect(ect2,bnds_plot);
nodes = [xplot yplot];
pdeplot( nodes', [], ect2c', 'xydata',hplot,'zdata',hplot, 'colormap', 'jet', 'mesh',mesh) ;
grid on
% xlabel('longitude')
% ylabel('latitude')
% zlabel('bathymetry (m)')
% title('linearly interpolated nodes')

% galveston trinity bay 
% caxis([.5 3.5])
% set(gca,'CameraPosition',[226105.1750978218 3500057.695551115 10.028497019894523])
% set(gca,'CameraTarget',[327500.0 3285000.0 2.0])
% set(gca,'CameraViewAngle',10.472173081209045)
% set(gcf,'PaperPosition',[.25 2.5 8 6])
% set(gcf,'Position',[1 1 1300 700])
% print('-r400','-dpng','~/galveston_linear_nodes.png') ;

%%
% Plot of rimls surface nodes
figure(5)
ect3 = delaunay(rimls_xplot,rimls_yplot);
ect3c = clean_ect(ect3,rimls_bnds_plot);
nodes = [rimls_xplot rimls_yplot];
pdeplot( nodes', [], ect3c', 'xydata',rimls_hplot,'zdata',rimls_hplot, 'colormap', 'jet', 'mesh',mesh) ;
grid on
% xlabel('longitude')
% ylabel('latitude')
% zlabel('bathymetry (m)')
% title('rimls surface nodes')

% galveston trinity bay
% caxis([.5 3.5])
% set(gca,'CameraPosition',[226105.1750978218 3500057.695551115 10.028497019894523])
% set(gca,'CameraTarget',[327500.0 3285000.0 2.0])
% set(gca,'CameraViewAngle',10.472173081209045)
% set(gcf,'PaperPosition',[.25 2.5 8 6])
% set(gcf,'Position',[1 1 1300 700])
% print('-r400','-dpng','~/galveston_rimls_nodes.png') ;

% % Plot of element normals
% figure(6)
% scale = 1000;
% nrm = nrm/scale;
% quiver3(cnds(:,1),cnds(:,2),cnds(:,3),nrm(:,1),nrm(:,2),nrm(:,3))



disp('Max absolute difference between rimls and  linear interpolation')
disp(max(abs(rimls_hpts-hpts)))
disp('Max relative difference between rimls and linear interpolation')
disp(max(abs(rimls_hpts-hpts))/max(hpts)*100)

disp('Max absolute difference between rimls and verticies')
disp(max(abs(rimls_vnds(:,3)-eval_hb)))
disp('Max relative difference between rimls and verticies')
disp(max(abs(rimls_vnds(:,3)-eval_hb))/max(eval_hb)*100)

disp('Max rimls value')
disp(max(rimls_hpts))
disp('Min rimls value')
disp(min(rimls_hpts))




% rmpath '/home/sbrus/Codes/SomeDGMaterials'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/devscript'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'
