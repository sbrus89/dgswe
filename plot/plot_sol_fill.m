restoredefaultpath
clear all
close all
clc

grd_direc = '~/Codes/dgswe/grids/';
sol_direc = '~/Codes/dgswe/work/';
% sol_direc = '~/data-drive/converge_quad/mesh1/P2/CTP2/';
% sol_direc = '/home/sbrus/data-drive/dgswe_converge_curve_bath/converge3/p3/ctp3/hbp3/';
% grd_name = 'inlet1.grd';
% grd_name = 'converge3.grd';
grd_name = 'converge1_dble.grd';
% grd_name = 'beaufort_hb+2.grd';
plot_folder = 'velplot';

% grd_direc = '~/Codes/dgswe/work/PE0000/';
% sol_direc = '~/Codes/dgswe/work/PE0000/';
% grd_name = 'fort.14';
% plot_folder = 'velplot';

% grd_direc = '~/Codes/dgswe/grids/';
% sol_direc = '~/data-drive/converge_quad/mesh1/P1/CTP2/';
% grd_name = 'converge_quad.grd';
% plot_folder = 'velplot';

% grd_direc = '~/data-drive/galveston/dgswe/tri_spline/p3/';
% sol_direc = '~/data-drive/galveston/dgswe/tri_spline/p3/';
% grd_name = 'galveston_tri_spline.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/tri/p3/';
% sol_direc = '~/data-drive/galveston/dgswe/tri/p3/';
% grd_name = 'galveston_tri.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/quad2/M2/p3/';
% sol_direc = '~/data-drive/galveston/dgswe/quad2/M2/p3/';
% grd_name = 'galveston2.grd';
% plot_folder = 'velplot_scale';
% 
% grd_direc = '~/data-drive/galveston/dgswe/quad2_spline/p1/';
% sol_direc = '~/data-drive/galveston/dgswe/quad2_spline/p1/';
% grd_name = 'galveston2_spline.grd';
% % grd_name = 'galveston2_plot.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/quad2_spline_channel/p3/';
% sol_direc = '~/data-drive/galveston/dgswe/quad2_spline_channel/p3/';
% grd_name = 'galveston2_spline_channel.grd';
% % grd_name = 'galveston2_plot.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/tri_ho-bath/p1/';
% sol_direc = '~/data-drive/galveston/dgswe/tri_ho-bath/p1/';
% grd_name = 'galveston_tri.grd';
% plot_folder = 'velplot_scale';

% grd_direc = '~/data-drive/galveston/dgswe/tri_ho-bath/p2/';
% sol_direc = '~/data-drive/galveston/dgswe/tri_ho-bath/p2/';
% grd_name = 'galveston_tri_spline.grd';
% plot_folder = 'velplot_scale';


ctp = 2;

nsnap = 49;

grid_on = 0;

% zoom_area = [3.2e5 3.4e5 3.24e6 3.26e6;
%              2.84e5 3.04e5 3.205e6 3.2325e6];
         
zoom_area = [];        

[nzoom,~] = size(zoom_area);

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,~] = size(EToV);
[nn,~] = size(VX);


% fid = fopen([grd_direc,'elem_nodes.d']);
% idx = zeros(ne,8);
% for i = 1:ne
%     line = fgetl(fid);
%     v = sscanf(line,' %d %d %f %f %f %f %f %f',8)';
%     idx(i,1) = v(1);
%     idx(i,2) = v(2);
%     idx(i,3:v(2)+2) = v(3:v(2)+2);
% end
% hbnodes = idx(:,3:end);


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

nvert = zeros(1,2*nel_type);
ndof = zeros(1,2*nel_type);
nnds = zeros(1,2*nel_type);
m2n = zeros(9,36,2*nel_type);
for et = 1:2*nel_type
    data  = fscanf(fid, '%d %d \n',2);
    ndof(et) = data(2);
    nvert(et) = data(1);    
    m2n(1:nvert(et),1:ndof(et),et) = fscanf(fid,' %g ',[ndof(et) nvert(et)])';
end
fclose(fid);

% ndof(3) = ndof(1);
% ndof(4) = ndof(2);

mndof = max(ndof);



fid_hb = fopen([sol_direc,'hb_modal.d']);
hbm = fscanf(fid_hb,'%g', [ne mndof])'; 
hb = zeros(9,ne,1);
for el = 1:ne
    
    n = nelnds(el);
    et = el_type(el);
    
    hb(1:n,el,1) = m2n(1:n,1:ndof(et+4),et+4)*hbm(1:ndof(et+4),el);
%     hb(1:n,el,1) = 10;
    
end
fclose(fid_hb);



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
Zv = zeros(9,ne,nsnap);
u = zeros(9,ne,nsnap);
v = zeros(9,ne,nsnap);

snap = 0;
while ~feof(fid_H) && snap < nsnap
    
    snap = snap + 1;
    
    th = fscanf(fid_H,' %g ', 1); % read in time
    Z = fscanf(fid_H,' %g ', [ne mndof])'; % read in H solution at time t
%     H = fscanf(fid_H,' %g ', [ne mndof])'; % read in H solution at time t    
    
    tqx = fscanf(fid_Qx,' %g ', 1); % read in time
    Qx = fscanf(fid_Qx,' %g ', [ne mndof])'; % read in Qx solution at time t
    
    t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
    Qy = fscanf(fid_Qy,' %g ', [ne mndof])'; % read in Qy solution at time t    
    
    
    for el = 1:ne
        
        n = nelnds(el);
        et = el_type(el);
        
        Zv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*Z(1:ndof(et),el);
%         Hv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*H(1:ndof(et),el);
        Qxv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*Qx(1:ndof(et),el);
        Qyv(1:n,el,snap) = m2n(1:n,1:ndof(et),et)*Qy(1:ndof(et),el);
        Hv(1:n,el,snap) = Zv(1:n,el,snap)+HB(EToV(el,1:n));
%         Hv(1:n,el,snap) = Zv(1:n,el,snap)+hbnodes(el,1:n)';
%         Zv(1:n,el,snap) = Hv(1:n,el,snap)-HB(EToV(el,1:n));
        
        u(1:n,el,snap) = Qxv(1:n,el,snap)./hb(1:n,el,1);
        v(1:n,el,snap) = Qyv(1:n,el,snap)./hb(1:n,el,1);
        
    end
    
end

% u = Qxv./Hv;
% v = Qyv./Hv;

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

figure

hold on

for el = 1:ne
    
    n = nelnds(el);
    et = el_type(el);
    
    if n == 9
        n = n-1;
    end
    
    if grid_on == 0
%         fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),HB(EToV(el,1:n)),'EdgeColor','none')
        fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),hb(1:n,el),'EdgeColor','none')
    else
%         fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),HB(EToV(el,1:n)))
        fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),hb(1:n,el))
    end
    
end
% 
% colorbar
% set(gcf,'PaperPositionMode','auto')
% axis equal
% hold off
% 
% print('-r350','-dpng',sprintf('%s/bathy',FramesFolder)) ;
% print('-depsc',sprintf('%s/bathy',FramesFolder)) ;
% 
% for zoom = 1:nzoom
%     axis(zoom_area(zoom,:))
%     print('-r350','-dpng',sprintf('%s/bathy_zoom%02d',FramesFolder,zoom)) ;
%     print('-depsc',sprintf('%s/bathy_zoom%02d',FramesFolder,zoom)) ;
% end
    
figure

for tsnap = 1:snap
    
%     figure('visible','off')
%     figure


    axis equal
    
    disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    [day,hr,minute,sec] = s2dhms(t(tsnap));
    
    hold on
%     velmin = min(min(vel(:,:,tsnap)));
%     velmax = max(max(vel(:,:,tsnap)));
%     caxis([velmin velmax])

    for el = 1:ne
        
        n = nelnds(el);
        et = el_type(el);
        
        if n == 9
            n = n-1;
        end
        
%         HB = 10 - 5*cos(2*pi/500*VX(EToV(el,1:n),2));
%         vel(1:n,el,tsnap) = sqrt((Qxv(1:n,el,tsnap)./HB).^2 + (Qyv(1:n,el,tsnap)./HB).^2);
        
        if grid_on == 0
            fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),vel(1:n,el,tsnap),'EdgeColor','none')            
        else            
            fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),vel(1:n,el,tsnap))
        end
        
%         quiver(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),u(1:n,el,tsnap),v(1:n,el,tsnap),'k')
    end
    
    hold off
    
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec),')'] ;
    ttext = ['Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec)] ;

    title(ttext)
    xlabel('x')
    ylabel('y')
    colorbar
    if scale == 1
        caxis([0 velscale(tsnap)])
    else
%         caxis([0 1])
%         caxis('auto')
    end
    axis image
    
    drawnow
    
    print('-r350','-dpng',sprintf('%s/vel%04d',FramesFolder,tsnap)) ;
%     print('-depsc',sprintf('%s/vel%04d',FramesFolder,tsnap)) ;
%     
%     for zoom = 1:nzoom
%         axis(zoom_area(zoom,:))
%         print('-r350','-dpng',sprintf('%s/vel%04d_zoom%02d',FramesFolder,tsnap,zoom)) ;
%         print('-depsc',sprintf('%s/vel%04d_zoom%02d',FramesFolder,tsnap,zoom)) ;
%     end
    
%         figure
%     
%         axis equal
%     
%         disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
%         [day,hr,minute,sec] = s2dhms(t(tsnap));
%     
%         hold on
%         
%         for el = 1:ne
%             
%             n = nelnds(el);
%             et = el_type(el);             
%             
%             if n == 9
%                 n = n-1;
%             end
%             
%             if grid_on == 0
%                 fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),Zv(1:n,el,tsnap),'EdgeColor','none')
% %                fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),hbnodes(el,1:n),'EdgeColor','none')
%             else
%                 fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),Zv(1:n,el,tsnap))
% %                fill(VX(EToV(el,1:n),1),VX(EToV(el,1:n),2),hbnodes(el,1:n))
%             end
%         end
%         
%         
%         hold off
%     
%         ttext = ['Surface solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(minute),', Second:  ',num2str(sec),')'] ;
%         title(ttext)
%         xlabel('x')
%         ylabel('y')
%         zmin = min(min(Zv(:,:,tsnap)));
%         zmax = max(max(Zv(:,:,tsnap)));
% %         caxis([Zmin Zmax])        
%         colorbar
%         axis image
%         
%         drawnow
%     
%     %     print('-r400','-dpng',sprintf('%s/vel%04d',FramesFolder,tsnap)) ;
%     %
%     %     for zoom = 1:nzoom
%     %         axis(zoom_area(zoom,:))
%     %         print('-r400','-dpng',sprintf('%s/vel%04d_zoom%02d',FramesFolder,tsnap,zoom)) ;
%     %     end
%     
        pause(.01)
end