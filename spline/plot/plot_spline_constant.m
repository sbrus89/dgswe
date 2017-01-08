% close all
clear all
clc

% out_direc = '../work/';
% out_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/';
out_direc = '/home/sbrus/data-drive/galveston_SL18/grid_dev/v24_cart/coarse/spline/';

file = fopen([out_direc,'spline.out']);

figure(1)

%     th = fscanf(fid_H,' %g ', 1); % read in time
%     Z = fscanf(fid_H,' %g ', [ne mndof])'; % read in H solution at time t

n_seg = fscanf(file,'%d',1);


for j = 1:n_seg
   n = fscanf(file, '%d',1)'; 
 

  
   
   coef_xy = fscanf(file, '%g',[10, n])';
   ax = coef_xy(:,1);
   bx = coef_xy(:,2);
   cx = coef_xy(:,3);
   dx = coef_xy(:,4);
   
   ay = coef_xy(:,5);
   by = coef_xy(:,6);
   cy = coef_xy(:,7);
   dy = coef_xy(:,8);   
   
   ti = coef_xy(:,9);
  
 
   t = linspace(ti(1),ti(n),10*n);
   
  
   %%%%%%%%%% Plot boundary polynomials %%%%%%%%%%

    x = zeros(size(t));
    y = zeros(size(t));
    xe = zeros(size(t));
    ye = zeros(size(t));
    for i = 1:n-1
        % Cubic Spline
        x = (t > ti(i) & t <= ti(i+1)).*(ax(i) + bx(i)*(t-ti(i)) + cx(i)*(t-ti(i)).^2 + dx(i)*(t-ti(i)).^3)+x;
        y = (t > ti(i) & t <= ti(i+1)).*(ay(i) + by(i)*(t-ti(i)) + cy(i)*(t-ti(i)).^2 + dy(i)*(t-ti(i)).^3)+y;
        
        % Linear element edge
        dt = ti(i+1)-ti(i);
        xe = (t > ti(i) & t <= ti(i+1)).*((t-ti(i))*ax(i+1)/dt - (t-ti(i)-dt)*ax(i)/dt) + xe;
        ye = (t > ti(i) & t <= ti(i+1)).*((t-ti(i))*ay(i+1)/dt - (t-ti(i)-dt)*ay(i)/dt) + ye;
    end
    xe(1) = ax(1);
    ye(1) = ay(1);
    xe(end) = ax(i+1);
    ye(end) = ay(i+1);
    
    x(1) = ax(1);
    y(1) = ay(1);
    x(end) = ax(i+1);
    y(end) = ay(i+1); 
    
%     figure(j)
    

    hold on
    plot(x,y,'-b','LineWidth',2)
    hold all
    plot(ax,ay,'o','MarkerSize',12,'MarkerFaceColor','b','Color','b')
%     plot(xe,ye,'k')

%     if j == 8
%     for i = 1:n
%        k = round(coef_xy(i,10));
%        txt = sprintf('%d',k) ;
%        text(ax(i),ay(i),txt,'FontSize',6) ; 
%     end
%     end
    

    hold off
end

fclose(file);




hold on

file = fopen([out_direc,'eval_nodes.out']);

n = textscan(file,'%d',1);
n_seg = n{1};

for j = 1:n_seg
%     figure(j)
    hold on
    
    n = textscan(file, '%d',1); 
    
    xy = textscan(file, '%f%f',n{1});
    

     plot(xy{:,1},xy{:,2},'x','MarkerSize',20,'Color','r')

end
fclose(file);




hold on

file = fopen([out_direc,'max_deform.out']);

n_seg = fscanf(file,'%d \n',1);

for j = 1:n_seg
%     figure(j)
    hold on
    
    nnds = textscan(file, '%d \n',1); 
    n = nnds{1};
    
    xy = fscanf(file, '%f',[7,n])';

    
     for i = 1:n
       plot([xy(i,1),xy(i,3)],[xy(i,2),xy(i,4)],'b','LineWidth',2)
     end
     

end
fclose(file);






% hold on
% 
% fid = fopen([out_direc,'fort.14_spline']) ;
% 
% agrid = fgetl(fid) ;
% disp(agrid) ;
% title = agrid ; 
% 
% N = fscanf(fid,'%g %g',2) ;
% 
% Val = zeros(N(2),4) ;
% Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;
% VX = Val(:,2:3) ;
% 
% plot(VX(:,1),VX(:,2),'mx')

axis equal
%axis([23000 27000 10000 12500])














