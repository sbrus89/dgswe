restoredefaultpath
clear all
close all
clc

% grd_direc = '~/Codes/dgswe/grids/';
% sol_direc = '~/Codes/dgswe/output/';
% % grd_name = 'inlet1_quad.grd';
% grd_name = 'converge_quad.grd';
% plot_folder = 'velplot_scale';

grd_direc = '~/data-drive/galveston/dgswe/quad2/M2/';
sol_direc = '~/data-drive/galveston/dgswe/quad2/M2/';
grd_name = 'galveston2.grd';
plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/quad2_spline/';
% sol_direc = '~/data-drive/galveston/dgswe/quad2_spline/';
% grd_name = 'galveston2_spline.grd';
% % grd_name = 'galveston2_plot.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/quad2_spline_channel/';
% sol_direc = '~/data-drive/galveston/dgswe/quad2_spline_channel/';
% grd_name = 'galveston2_spline_channel.grd';
% % grd_name = 'galveston2_plot.grd';
% plot_folder = 'velplot_scale';

ctp = 2;

nsnap = 200;

grid_on = 1;

zoom_area = [3.2e5 3.4e5 3.24e6 3.26e6;
             2.84e5 3.04e5 3.205e6 3.2325e6];
         
% zoom_area = [];        

[nzoom,~] = size(zoom_area);

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,~] = size(EToV);
[nn,~] = size(VX);

FramesFolder = strcat(sol_direc,plot_folder) ;
if ( exist(FramesFolder,'dir') == 0 ) 
    mkdir(FramesFolder) ; 
end

el_type = zeros(ne,1);
for el = 1:ne
   if nelnds(el) == 3
       el_type(el) = 1;
   elseif nelnds(el) == 4
       el_type(el) = 2;
   elseif nelnds(el) == (ctp+1)*(ctp+2)/2
       el_type(el) = 3;
   elseif nelnds(el) == (ctp+1)^2
       el_type(el) = 4;
   end
       
end

fid = fopen([sol_direc,'modal2nodal.d']);

nel_type = 4;

nvert = zeros(1,nel_type);
ndof = zeros(1,nel_type);
m2n = zeros(4,36,nel_type);
for et = 1:nel_type
    data  = fscanf(fid, '%d %d \n',2);
    ndof(et) = data(2);
    nvert(et) = data(1);
    m2n(1:nvert(et),1:ndof(et),et) = fscanf(fid,' %g ',[ndof(et) nvert(et)])';
end
fclose(fid);

mndof = max(ndof);

fid_H = fopen([sol_direc,'solution_H.d']);
fid_Qx = fopen([sol_direc,'solution_Qx.d']);
fid_Qy = fopen([sol_direc,'solution_Qy.d']);

line = fgetl(fid_H);
line = fgetl(fid_Qx);
line = fgetl(fid_Qy);

% Hv = zeros(4,ne,nsnap);
% Qxv = zeros(4,ne,nsnap);
% Qyv = zeros(4,ne,nsnap);
% zv = zeros(4,ne,nsnap);
Hv = zeros(9,ne,nsnap);
Qxv = zeros(9,ne,nsnap);
Qyv = zeros(9,ne,nsnap);
zv = zeros(9,ne,nsnap);

snap = 0;
while ~feof(fid_H) && snap < nsnap
    
    snap = snap + 1;
    
    th = fscanf(fid_H,' %g ', 1); % read in time
    H = fscanf(fid_H,' %g ', [ne mndof])'; % read in H solution at time t
    
    tqx = fscanf(fid_Qx,' %g ', 1); % read in time
    Qx = fscanf(fid_Qx,' %g ', [ne mndof])'; % read in Qx solution at time t
    
    t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
    Qy = fscanf(fid_Qy,' %g ', [ne mndof])'; % read in Qy solution at time t
    
    %     for el = 1:ne
    %         if nelnds(el) == 3
    %             Hv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*H(1:ndof(1),el);
    %             Qxv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*Qx(1:ndof(1),el);
    %             Qyv(1:3,el,snap) = m2n(1:3,1:ndof(1),1)*Qy(1:ndof(1),el);
    %             zv(1:3,el,snap) = Hv(1:3,el,snap)-HB(EToV(el,1:3));
    %         elseif nelnds(el) == 4
    %             Hv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*H(1:ndof(2),el);
    %             Qxv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*Qx(1:ndof(2),el);
    %             Qyv(1:4,el,snap) = m2n(1:4,1:ndof(2),2)*Qy(1:ndof(2),el);
    %             zv(1:4,el,snap) = Hv(1:4,el,snap)-HB(EToV(el,1:4));
    %         end
    %     end
    
    for el = 1:ne
        
        n = nelnds(el);
        et = el_type(el);
        
        Hv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*H(1:ndof(et),el);
        Qxv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*Qx(1:ndof(et),el);
        Qyv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*Qy(1:ndof(et),el);
        zv(1:n,el,snap) = Hv(1:n,el,snap)-HB(EToV(el,1:n));
        
    end
end

u = Qxv./Hv;
v = Qyv./Hv;
vel = sqrt(u.^2 + v.^2);

Hmax = max(max(max(Hv)));
Qxmax = max(max(max(Qxv)));
Qymax = max(max(max(Qyv)));
velmax = max(max(max(vel)));

Hmin = min(min(min(Hv)));
Qxmin = min(min(min(Qxv)));
Qymin = min(min(min(Qyv)));
velmin = min(min(min(vel)));

scale = 0;
if exist([sol_direc,'velscale.mat'])
    load([sol_direc,'velscale.mat'])
    disp('velscale.mat was found and loaded')
    scale = 1;
end


for tsnap = snap %1:snap
    
    figure
    
    axis equal
    
    disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    [day,hr,minute,sec] = s2dhms(t(tsnap));
    
    hold on
    velmin = min(min(vel(:,:,tsnap)));
    velmax = max(max(vel(:,:,tsnap)));
    caxis([velmin velmax])
    %     for el = 1:ne
    %
    %         if grid_on == 0
    %             if nelnds(el) == 3
    %                 fill(VX(EToV(el,1:3),1),VX(EToV(el,1:3),2),vel(1:3,el,tsnap),'EdgeColor','none')
    %             elseif nelnds(el) == 4
    %                 fill(VX(EToV(el,1:4),1),VX(EToV(el,1:4),2),vel(1:4,el,tsnap),'EdgeColor','none')
    %             end
    %         else
    %             if nelnds(el) == 3
    %                 fill(VX(EToV(el,1:3),1),VX(EToV(el,1:3),2),vel(1:3,el,tsnap))
    %             elseif nelnds(el) == 4
    %                 fill(VX(EToV(el,1:4),1),VX(EToV(el,1:4),2),vel(1:4,el,tsnap))
    %             end
    %         end
    %     end
    
    for el = 1:ne
        
        n = nelnds(el);
        et = el_type(el);
        
        if n == 9
            n = n-1;
        end
        
        if grid_on == 0
            fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),vel(1:n,el,tsnap),'EdgeColor','none')            
        else            
            fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),vel(1:n,el,tsnap))
        end
        
%         quiver(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),u(1:n,el,tsnap),v(1:n,el,tsnap),'k')
    end
    
    hold off
    
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    colorbar
    if scale == 1
        caxis([0 velscale(tsnap)])
    else
        caxis([0 1])
    end
    axis image
    
    pause(.01)
    
    print('-r350','-dpng',sprintf('%s/vel%04d',FramesFolder,tsnap)) ;
    
    for zoom = 1:nzoom
        axis(zoom_area(zoom,:))
        print('-r350','-dpng',sprintf('%s/vel%04d_zoom%02d',FramesFolder,tsnap,zoom)) ;
    end
    
    %     figure
    %
    %     axis equal
    %
    %     disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    %     [day,hr,minute,sec] = s2dhms(t(tsnap));
    %
    %     hold on
    %     zmin = min(min(zv(:,:,tsnap)));
    %     zmax = max(max(zv(:,:,tsnap)));
    %     caxis([zmin zmax])
    %     for el = 1:ne
    %
    %         if grid_on == 0
    %             if nelnds(el) == 3
    %                 fill(VX(EToV(el,1:3),1),VX(EToV(el,1:3),2),zv(1:3,el,tsnap),'EdgeColor','none')
    %             elseif nelnds(el) == 4
    %                 fill(VX(EToV(el,1:4),1),VX(EToV(el,1:4),2),zv(1:4,el,tsnap),'EdgeColor','none')
    %             end
    %         else
    %             if nelnds(el) == 3
    %                 fill(VX(EToV(el,1:3),1),VX(EToV(el,1:3),2),zv(1:3,el,tsnap))
    %             elseif nelnds(el) == 4
    %                 fill(VX(EToV(el,1:4),1),VX(EToV(el,1:4),2),zv(1:4,el,tsnap))
    %             end
    %         end
    %     end
    %     hold off
    %
    %     ttext = ['Surface solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec),')'] ;
    %     title(ttext)
    %     xlabel('x')
    %     ylabel('y')
    %     colorbar
    %     axis image
    
    %     print('-r400','-dpng',sprintf('%s/vel%04d',FramesFolder,tsnap)) ;
    %
    %     for zoom = 1:nzoom
    %         axis(zoom_area(zoom,:))
    %         print('-r400','-dpng',sprintf('%s/vel%04d_zoom%02d',FramesFolder,tsnap,zoom)) ;
    %     end
    
    %     pause(.01)
end