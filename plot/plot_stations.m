clc
close all
clear all

adc_sta_direc = {
%                 '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/adcirc/ESL0_orig/';
                '/home/sbrus/data-drive/galveston_spline_flux_x10/galveston_tri/adcirc/ESL0/';
                '/home/sbrus/data-drive/galveston_spline_flux_x10/galveston_tri_x64/adcirc/ESL0/';
                };         

dg_sta_direc = {
               '/home/sbrus/data-drive/galveston_spline_flux_x10/galveston_tri_x64/p3/ctp3/hbp3/';
               '/home/sbrus/data-drive/galveston_spline_flux_x10/galveston_tri/p3/ctp3/hbp3/';
%                '/home/sbrus/data-drive/galveston_spline_flux_x10/galveston_tri/p3/ctp1/hbp1/'               
               };         
         
        
% grid_name = '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p1/galveston2_plot.grd';

grid_name = '/home/sbrus/data-drive/galveston_spline_flux_x10/grids/galveston_tri.grd';
stations_file = '/home/sbrus/data-drive/galveston_spline_flux_x10/grids/stations.d';
save_dir = '/home/sbrus/data-drive/galveston_spline_flux_x10/';         
         
cmap = colormap(lines);          

labels = {'adcirc';
          'adcirc x64';
          'DG x64 p3 tri h.o.';
          %'DG p2 tri';
          'DG p3 tri h.o.'
          %'DG p3 tri lin'
          %'DG p1 quad';
          %'DG p2 quad'
          %'DG p3 quad'
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

% nsnap_plot = 552;
% ts_sta = 100;

nsnap_read = 576;
nsta = 601;

% xs_snap_vec = 2:576;
% ts_sta_vec = 1:481;

xs_snap_vec = 2:576;
ts_sta_vec = 1:481;


fig_onoff = 'off';
fig_size = [0 0 11 4.5];
fig_fmt = 'png';
res = '300';


nruns_dg = length(dg_sta_direc);
nruns_adc = length(adc_sta_direc);


fid = fopen(stations_file);
nsta = fscanf(fid,' %g ', 1);
xysta = fscanf(fid,' %g ', [2 nsta])';

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
        i = i+1;
    end
    set(gcf, 'units', 'inches', 'pos', fig_size);
    set(gcf,'PaperPositionMode','auto') ;    
    legend(labels,'Location','Northwest','FontSize',10) 
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
        i = i+1;        
    end
    set(gcf, 'units', 'inches', 'pos', fig_size)
    set(gcf,'PaperPositionMode','auto') ;    
    legend(labels,'Location','Northwest','FontSize',10)
    title(['Center channel stations, t = ',num2str(t(snap)),' s'])    
    xlabel('station')  
    ylabel('velocity (m/s)')
    print(['-d',fig_fmt],['-r',res],[save_dir,'vel_xs_snap',num2str(snap),'.',fig_fmt])   
    hold off
    
    j = j+1;
end


% Time series plots
j = 1;
for sta = ts_sta_vec
    disp(['  station: ',num2str(sta), ', # ',num2str(j),'/',num2str(length(ts_sta_vec))])
    
    figure(3)
    set(gcf,'visible',fig_onoff) 
    subplot(1,4,[1 3])    
    i = 1;
    for run = 1:nruns_adc
        plot(time,z_adc(sta,:,run),'Color',colors(i,:))
        hold on
        i = i+1;
    end
    for run = 1:nruns_dg
        plot(t,z_dg(sta,:,run),'Color',colors(i,:))
        i = i+1;        
    end
    set(gcf, 'units', 'inches', 'pos', fig_size);
    set(gcf,'PaperPositionMode','auto') ;    
    legend(labels,'Location','Northwest','FontSize',10);
    title(['Time series at station: ',num2str(sta)])      
    xlabel('t (s)')
    ylabel('zeta (m)')
    hold off
    subplot(1,4,4);
    sta_dot = plot(xysta(sta,1),xysta(sta,2),'go','MarkerFaceColor','g','MarkerSize',5);    
    print(['-d',fig_fmt],['-r',res],[save_dir,'zeta_ts_sta_',num2str(sta),'.',fig_fmt]) 
    delete(sta_dot)
    
    figure(4)
    set(gcf,'visible',fig_onoff)
    subplot(1,4,[1 3]) 
    i = 1;
    for run = 1:nruns_adc
        plot(time,vel_adc(sta,:,run),'Color',colors(i,:))
        hold on
        i = i+1;
    end
    for run = 1:nruns_dg
        plot(t,vel_dg(sta,:,run),'Color',colors(i,:))
        i = i+1;        
    end
    set(gcf, 'units', 'inches', 'pos', fig_size)
    set(gcf,'PaperPositionMode','auto') ;
    legend(labels,'Location','Northwest','FontSize',10) 
    title(['Time series at station: ',num2str(sta)])     
    xlabel('t (s)')
    ylabel('velocity (m/s)')    
    hold off
    subplot(1,4,4);
    sta_dot = plot(xysta(sta,1),xysta(sta,2),'go','MarkerFaceColor','g','MarkerSize',5);     
    print(['-d',fig_fmt],['-r',res],[save_dir,'vel_ts_sta_',num2str(sta),'.',fig_fmt])
    delete(sta_dot)
    
    
    j = j+1;
end




% % quiver(xyH(sta,1),xyH(sta,2),u(sta,snap),v(sta,snap),'b')
