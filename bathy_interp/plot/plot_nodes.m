clear all
close all
clc

% path(path,'/home/sbrus/Codes/SomeDGMaterials')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/devscript')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
% path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

out_dir = '../output/';

finp = fopen('../work/bathy.inp');

n = 0;
while n < 8 && ~feof(finp)
   temp = strtrim(fgetl(finp));
     
   if isempty(temp) || temp(1) == '!'
       
   else
      n = n + 1;
      temp2 = strtrim(strsplit(temp));
      names{n} = temp2{1};
   end
end

base_name = names{1};
eval_name = names{6};


mesh = 'on';

% entire area
xbox = [0 0];
ybox = [0 0];

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


[base_ect,base_xy,base_hb,nelnds,opedat,boudat,title,base_bvnds] = readfort14(base_name);
if strcmp(base_name,eval_name)
    eval_ect = base_ect;
    eval_xy = base_xy;
    eval_hb = base_hb;
    eval_bvnds = base_bvnds;
else
    [eval_ect,eval_xy,eval_hb,nelnds,opedat,boudat,title,eval_bvnds] = readfort14(eval_name);
end


[base_ne,~] = size(base_ect);
[eval_ne,~] = size(eval_ect);


fid = fopen([out_dir,'boundary_nodes.d']);
N = fscanf(fid,'%g',1) ;
bnd_flag = fscanf(fid,'%g \n', [1 N(1)])' ;
fclose(fid);

fid = fopen([out_dir,'interp_nodes.d']);
N = fscanf(fid,'%g',1) ;
nodes = fscanf(fid,'%g %g %g \n', [3 N(1)])' ;
fclose(fid);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
hnodes = nodes(:,3);

if all(xbox) == 0 && all(ybox) == 0
    xbox = [min(xnodes) max(xnodes)];
    ybox = [min(ynodes) max(ynodes)];
end

in = inpolygon(xnodes,ynodes,xbox,ybox);

xplot = xnodes(in);
yplot = ynodes(in);
bnds_plot = bnd_flag(in);


% Plot of full grid batymetry and zoom box
if eval_ne > base_ne
    figure(1)
    trimesh(eval_ect,eval_xy(:,1),eval_xy(:,2),0*eval_xy(:,1),'EdgeColor',[1 0 0],'LineStyle','-.','LineWidth',2)
    hold on
    trimesh(base_ect,base_xy(:,1),base_xy(:,2),0*base_xy(:,1),'EdgeColor',[0 0 1],'LineWidth',2)
else
    trimesh(base_ect,base_xy(:,1),base_xy(:,2),0*base_xy(:,1),'EdgeColor',[0 0 1],'LineWidth',2)
    hold on
    trimesh(eval_ect,eval_xy(:,1),eval_xy(:,2),0*eval_xy(:,1),'EdgeColor',[1 0 0],'LineStyle','-.','LineWidth',2)
end

plot(xnodes,ynodes,'r.','MarkerSize',20)
plot([xbox(1) xbox(2)],[ybox(1) ybox(1)],'r','LineWidth',5)
plot([xbox(1) xbox(2)],[ybox(2) ybox(2) ],'r','LineWidth',5)
plot([xbox(1) xbox(1)],[ybox(1) ybox(2) ],'r','LineWidth',5)
plot([xbox(2) xbox(2)],[ybox(1) ybox(2) ],'r','LineWidth',5)
xlabel('longitude')
ylabel('latitude')
view(0,90)
axis image
grid on


% Plot of interpolated nodes
figure(4)
ect = delaunay(xplot,yplot);
ectc = clean_ect(ect,bnds_plot);
nodes = [xplot yplot];
pdeplot( nodes', [], ectc', 'xydata',hnodes,'zdata',hnodes, 'colormap', 'jet', 'mesh',mesh) ;
grid on
xlabel('longitude')
ylabel('latitude')
zlabel('bathymetry (m)')






% rmpath '/home/sbrus/Codes/SomeDGMaterials'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/devscript'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
% rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'
