clc
close all
clear all

adc_sta_direc = {
%                 '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0_orig/';
%                 '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0/';
%                 '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x4/adcirc/ESL0/';
%                 '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x16/adcirc/ESL0/';                
%                 '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/';
% 
%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri/adcirc/';  
%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri_x64/adcirc/';      

% %                   '/home/sbrus/data-drive/galveston_SL18_tides/galveston_SL18_cart/esl1/adcirc_default/';
% %                   '/home/sbrus/data-drive/galveston_SL18_tides/galveston_SL18_cart/esl1/adcirc_513111/';     
%                   '/home/sbrus/data-drive/galveston_SL18_tides/galveston_SL18_cart/esl1/adcirc/';                  
% %                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse/adcirc/';
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl1/adcirc/';
                };         

dg_sta_direc = {
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/';
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x16/p3/ctp3/hbp3/';
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x4/p3/ctp3/hbp3/';              
               '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/';
               '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/rk22/';               
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp1/hbp1/'                      
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p2/ctp2/hbp2/';
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p1/ctp2/hbp1/';

%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x64/p3/ctp3/hbp3/';
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x16/p3/ctp3/hbp3/';
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad_x4/p3/ctp3/hbp3/';              
%                '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_quad/p3/ctp3/hbp3/';

%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri/p1/ctp1/hbp1/'; 
%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri/p1/ctp2/hbp1/';
%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri/p3/ctp3/hbp3/'
%                '/home/sbrus/data-drive/galveston_spline_surge/galveston_tri/p2/ctp2/hbp2/'

%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse/esl0/p1/ctp2/hbp1/';    
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl0/p3/ctp3/hbp3/';    
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl.5/p3/ctp3/hbp3/rk45/';
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl.5/p3/ctp1/hbp1/rk45/';     
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl.5/p3/ctp3/hbp3/';
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl.5/p3/ctp3/hbp3/rk22_str_tol0/';                  
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl.5/p3/ctp1/hbp1/';                     
%                   '/home/sbrus/data-drive/galveston_SL18_tides/coarse_x2/esl1/p3/ctp3/hbp3/'; 

%                   '/home/sbrus/data-drive/galveston_SL18_surge/coarse_x2/esl.5/p3/ctp3/hbp3/';
%                   '/home/sbrus/data-drive/galveston_SL18_surge/coarse_x2/esl.5/p3/ctp1/hbp1/';
               };         
         
        
grid_name = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_tri.grd';
% grid_name = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_quad.grd';
stations_file = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/stations.d';
save_dir = '/home/sbrus/data-drive/galveston_spline_flux_fix/';         

% grid_name = '/home/sbrus/data-drive/galveston_spline_surge/grids/galveston_tri.grd';
% stations_file = '/home/sbrus/data-drive/galveston_spline_surge/grids/stations.d';
% save_dir = '/home/sbrus/data-drive/galveston_spline_surge/';   

% grid_name = '/home/sbrus/data-drive/galveston_SL18_tides/grids/coarse_x2.grd';
% stations_file = '/home/sbrus/data-drive/galveston_SL18_tides/grids/stations_avg.d';
% save_dir = '/home/sbrus/data-drive/galveston_SL18_tides/';   

% grid_name = '/home/sbrus/data-drive/galveston_SL18_surge/grids/coarse_x2.grd';
% stations_file = '/home/sbrus/data-drive/galveston_SL18_surge/grids/stations_avg.d';
% save_dir = '/home/sbrus/data-drive/galveston_SL18_surge/';   
         
cmap = colormap(lines);          

labels = {
%           'adcirc';
%           'adcirc x4';
%           'adcirc x16';          
%           'adcirc x64';
%             'DG p1 tri lin';
%             'DG p1 tri h.o.'            
%           'DG x64 p3 tri h.o.';
%           'DG x16 p3 tri h.o.';
%           'DG x4 p3 tri h.o.';
%             'DG p3 tri h.o.';
%             'DG p2 tri h.o.';
%           'DG p3 tri lin'            

%           'DG x64 p3 quad h.o.'
%           'DG x16 p3 quad h.o.'
%           'DG x4 p3 quad h.o.'
%           'DG p3 quad h.o.'

% %             'adcirc fine default';
% %             'adcirc fine 513111'
%             'adcirc fine';            
% %             'adcirc med';
%             'adcirc coarse'
            
% %             'dg med p=1'
%             'dg coarse p=3 esl=0'
%             'dg coarse p=3 strtol=80'
%             'dg coarse p=3 h.o.'
%             'dg coarse p=3 lin'            
% %             'dg coarse p=3 esl=1'
               'dg h.o. rk45';
%                'dg lin rk45';
               'dg h.o. rk22';
%                'dg lin rk 22'
          }; 
colors = [cmap(1,:)
          cmap(2,:)
          cmap(3,:)
          cmap(4,:)
          cmap(5,:)];      
     



sta_vec = 1:481; % center of channel
% sta_vec = 482:521;
% sta_vec = 522:561;
% sta_vec = 542:601;

% sta_vec = 1:1763;

% nsnap_plot = 552;
% ts_sta = 100;

nsnap_read = 576;

xs_flag = 1;
ts_flag = 1;
sta_loc_flag = 0;

% xs_snap_vec = 2:576;
% ts_sta_vec = 1:481;

xs_snap_vec = 2:10:nsnap_read;
ts_sta_vec = 1:20:sta_vec(end);


fig_onoff = 'off';
fig_size = [0 0 11 4.5];
fig_fmt = 'png';
res = '300';


nruns_dg = length(dg_sta_direc);
nruns_adc = length(adc_sta_direc);


fid = fopen(stations_file);
nsta = fscanf(fid,' %g ', 1);
xysta = fscanf(fid,' %g ', [2 nsta])';

if sta_loc_flag == 1
    figure(3)
    set(gcf,'visible',fig_onoff)
    subplot(1,4,4)
    hold on
    plot_grid_function(grid_name);
    plot(xysta(:,1),xysta(:,2),'ro','MarkerFaceColor','r','MarkerSize',2)
    axis off
    
    figure(4)
    set(gcf,'visible',fig_onoff)
    subplot(1,4,4)
    hold on
    plot_grid_function(grid_name);
    plot(xysta(:,1),xysta(:,2),'ro','MarkerFaceColor','r','MarkerSize',2)
    axis off
end

sta_dist = zeros(nsta);
for sta = 1:nsta-1
   sta_dist(sta+1) = sta_dist(sta) + sqrt((xysta(sta+1,1)-xysta(sta,1))^2 + (xysta(sta+1,2)-xysta(sta,2))^2);
end



z_adc = zeros(nsta,nsnap_read,nruns_adc);
vel_adc = zeros(nsta,nsnap_read,nruns_adc);
for run = 1:nruns_adc
    [z,vel,time] = read_adcirc_stations(adc_sta_direc{run});
    
    z_adc(:,:,run) = z(:,:);
    vel_adc(:,:,run) = vel(:,:);
    t = time;
end


z_dg = zeros(nsta,nsnap_read,nruns_dg);
vel_dg = zeros(nsta,nsnap_read,nruns_dg);
for run = 1:nruns_dg    
    [Z,t] = read_solution(dg_sta_direc{run},'Z.sta',nsnap_read);
    [Qx,~] = read_solution(dg_sta_direc{run},'Qx.sta',nsnap_read);
    [Qy,~] = read_solution(dg_sta_direc{run},'Qy.sta',nsnap_read);
    [hb,~] = read_solution(dg_sta_direc{run},'hb.sta',nsnap_read);
    
    Z = squeeze(Z);
    Qx = squeeze(Qx);
    Qy = squeeze(Qy);
    hb = squeeze(hb);
    
    H = zeros(size(Z));
    
    for snap = 1:nsnap_read
        H(:,snap) = Z(:,snap) + hb;
    end
    
    Vel = sqrt((Qx./H).^2 + (Qy./H).^2);
    
    z_dg(:,:,run) = Z(:,:);
    vel_dg(:,:,run) = Vel(:,:);
end


% Cross section plots    
if xs_flag == 1
    j = 1;
    for snap = xs_snap_vec
        disp(['  xsection: ',num2str(snap), ', # ',num2str(j),'/',num2str(length(xs_snap_vec))])
        
        
        figure(1)
        set(gcf,'visible',fig_onoff)
        i = 1;
        for run = 1:nruns_adc
            plot(z_adc(sta_vec,snap-1,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        for run = 1:nruns_dg
            plot(z_dg(sta_vec,snap,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        set(gcf, 'units', 'inches', 'pos', fig_size);
        set(gcf,'PaperPositionMode','auto') ;
        legend(labels,'Location','Northeast','FontSize',10)
        title(['Center channel stations, t = ',num2str(t(snap)),' s'])
        xlabel('station')
        ylabel('zeta (m)')
        print(['-d',fig_fmt],['-r',res],[save_dir,'zeta_xs_snap_',num2str(snap),'.',fig_fmt])
        hold off
        
        figure(2)
        set(gcf,'visible',fig_onoff)
        i = 1;
        for run = 1:nruns_adc
            plot(vel_adc(sta_vec,snap-1,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        for run = 1:nruns_dg
            plot(vel_dg(sta_vec,snap,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        set(gcf, 'units', 'inches', 'pos', fig_size)
        set(gcf,'PaperPositionMode','auto') ;
        legend(labels,'Location','Northeast','FontSize',10)
        title(['Center channel stations, t = ',num2str(t(snap)),' s'])
        xlabel('station')
        ylabel('velocity (m/s)')
        print(['-d',fig_fmt],['-r',res],[save_dir,'vel_xs_snap',num2str(snap),'.',fig_fmt])
        hold off
        
        j = j+1;
    end
end

% Time series plots
if ts_flag == 1
    j = 1;
    for sta = ts_sta_vec
        disp(['  station: ',num2str(sta), ', # ',num2str(j),'/',num2str(length(ts_sta_vec))])
        
        figure(3)
        set(gcf,'visible',fig_onoff)
        if sta_loc_flag == 1
            subplot(1,4,[1 3])
        end
        i = 1;
        for run = 1:nruns_adc
            plot(time,z_adc(sta,:,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        for run = 1:nruns_dg
            plot(t,z_dg(sta,:,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        set(gcf, 'units', 'inches', 'pos', fig_size);
        set(gcf,'PaperPositionMode','auto') ;
        legend(labels,'Location','Northwest','FontSize',10);
        title(['Time series at station: ',num2str(sta)])
        xlabel('t (s)')
        ylabel('zeta (m)')
        hold off
        if sta_loc_flag == 1
            subplot(1,4,4);
            sta_dot = plot(xysta(sta,1),xysta(sta,2),'go','MarkerFaceColor','g','MarkerSize',5);
        end
        
            print(['-d',fig_fmt],['-r',res],[save_dir,'zeta_ts_sta_',num2str(sta),'.',fig_fmt])
            
        if sta_loc_flag == 1
            delete(sta_dot)
        end
        
        figure(4)
        set(gcf,'visible',fig_onoff)
        if sta_loc_flag == 1            
            subplot(1,4,[1 3])
        end
        i = 1;
        for run = 1:nruns_adc
            plot(time,vel_adc(sta,:,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        for run = 1:nruns_dg
            plot(t,vel_dg(sta,:,run),'Color',colors(i,:))
            hold on
            i = i+1;
        end
        set(gcf, 'units', 'inches', 'pos', fig_size)
        set(gcf,'PaperPositionMode','auto') ;
        legend(labels,'Location','Northwest','FontSize',10)
        title(['Time series at station: ',num2str(sta)])
        xlabel('t (s)')
        ylabel('velocity (m/s)')
        hold off
        if sta_loc_flag == 1
            subplot(1,4,4);
            sta_dot = plot(xysta(sta,1),xysta(sta,2),'go','MarkerFaceColor','g','MarkerSize',5);
        end
        
            print(['-d',fig_fmt],['-r',res],[save_dir,'vel_ts_sta_',num2str(sta),'.',fig_fmt])
            
        if sta_loc_flag == 1
            delete(sta_dot)
        end
        
        
        j = j+1;
    end
end



% % quiver(xyH(sta,1),xyH(sta,2),u(sta,snap),v(sta,snap),'b')
