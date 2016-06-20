clc
close all
clear all

         
sta_direc = {%'/home/sbrus/data-drive/galveston_spline/tri/p3/ctp3/hbp3/';
             %'/home/sbrus/data-drive/galveston_spline/tri/p3/ctp1/hbp1/';
%              %'/home/sbrus/data-drive/galveston_spline/quad/p3/ctp3/hbp3/';
%              '/home/sbrus/data-drive/galveston_spline/tri_x64/p1/ctp1/hbp1/';
%              '/home/sbrus/data-drive/galveston_spline/tri_x64/p1/ctp2/hbp1/';
%              '/home/sbrus/data-drive/galveston_spline/tri_x64/p1/ctp3/hbp1/';             
               '/home/sbrus/Codes/dgswe/work/';
             };         
         
adc_direc1 = '/home/sbrus/data-drive/galveston_spline/adcirc/tri/ESL0/';
adc_direc2 = '/home/sbrus/data-drive/galveston_spline/adcirc/tri_x64/ESL0/';
        
% grid_name = '/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/p1/galveston2_plot.grd';
grid_name = '/home/sbrus/data-drive/galveston_spline/grids/galveston_tri.grd';
         
         
cmap = colormap(lines);          

labels = {'adcirc';
          'adcircx64';
          %'DG p1 tri';
          %'DG p2 tri';
          'DG p3 tri'
          %'DG p1 quad';
          %'DG p2 quad'
          'DG p3 quad'}; 
colors = [cmap(1,:)
          cmap(2,:)
          cmap(3,:)
          cmap(4,:)
          cmap(5,:)];      
     



sta = 1:481; % center of channel
% sta = 482:521;
% sta = 522:561;
% sta = 542:601;

nsnap = 553;
ts_sta = 50;







nruns = length(sta_direc);

figure(1)
hold all
figure(2)
hold all
figure(3)
hold all
figure(4)
hold all

plot_adcirc_stations(adc_direc1,nsnap,sta,ts_sta,colors(1,:))
plot_adcirc_stations(adc_direc2,nsnap,sta,ts_sta,colors(2,:))






for run = 1:nruns
    

[Z,t] = read_solution(sta_direc{run},'Z.sta',nsnap);
[Qx,~] = read_solution(sta_direc{run},'Qx.sta',nsnap);
[Qy,~] = read_solution(sta_direc{run},'Qy.sta',nsnap);
[hb,~] = read_solution(sta_direc{run},'hb.sta',nsnap);

Z = squeeze(Z);
Qx = squeeze(Qx);
Qy = squeeze(Qy);
hb = squeeze(hb);

H = zeros(size(Z));

for snap = 1:nsnap    
    H(:,snap) = Z(:,snap) + hb;
end


u = Qx./H;
v = Qy./H;


% nsnap = snap;

for snap = nsnap
    figure(1)
    plot(Z(sta,snap),'Color',colors(run+2,:));
    figure(2)
    plot(sqrt(u(sta,snap).^2 + v(sta,snap).^2),'Color',colors(run+2,:))

end

figure(3)
plot(t,Z(ts_sta,:),'Color',colors(run+2,:));
figure(4)
plot(t,sqrt(u(ts_sta,:).^2 + v(ts_sta,:).^2),'Color',colors(run+2,:))

end





figure
hold on
plot_grid_function(grid_name);
% fid = fopen([sta_direc{run},'stations.d']);
fid = fopen('/home/sbrus/data-drive/galveston_spline/grids/stations.d');
nsta = fscanf(fid,' %g ', 1);
xysta = fscanf(fid,' %g ', [2 nsta])';
plot(xysta(:,1),xysta(:,2),'ro','MarkerFaceColor','r','MarkerSize',2)
plot(xysta(ts_sta,1),xysta(ts_sta,2),'bo','MarkerFaceColor','b','MarkerSize',7)


% figure(1)
% % legend(labels,'Location','Southwest','FontSize',10)
% xlabel('station')
% ylabel('surface elevation (m)')
% set(gcf, 'units', 'inches', 'pos', [0 0 3.5 2.8])
% set(gcf,'PaperPositionMode','auto') ;
% print('-dpdf',['~/Presentations/proposal/figures/galveston_zeta_',num2str(nsnap),'.pdf'])
% figure(2)
% % legend(labels,'Location','Northwest','FontSize',10)
% xlabel('station')
% ylabel('velocity (m/s)')
% set(gcf, 'units', 'inches', 'pos', [0 0 3.5 2.8])
% set(gcf,'PaperPositionMode','auto') ;
% print('-dpdf',['~/Presentations/proposal/figures/galveston_vel_',num2str(nsnap),'.pdf'])


% % quiver(xyH(sta,1),xyH(sta,2),u(sta,snap),v(sta,snap),'b')
