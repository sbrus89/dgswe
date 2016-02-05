close all
clear all
clc

t = 0:.00001:1;

file = fopen('work/spline.out');

figure

n = textscan(file,'%d',1);
n_seg = n{1};

for j = 1:n_seg
   n_dt = textscan(file, '%d%f',1); 
   n_seg_nodes(j) = n_dt{1,1};
   dt(j) = n_dt{1,2};
  
   
   coef_x = textscan(file, '%f%f%f%f',n_seg_nodes(j));
   ax = coef_x{:,1};
   bx = coef_x{:,2};
   cx = coef_x{:,3};
   dx = coef_x{:,4};
   
   coef_y = textscan(file, '%f%f%f%f',n_seg_nodes(j));
   ay = coef_y{:,1};
   by = coef_y{:,2};
   cy = coef_y{:,3};
   dy = coef_y{:,4};
   
    
   %%%%%%%%%% Plot boundary polynomials %%%%%%%%%%
    ti = 0;
    x = zeros(size(t));
    y = zeros(size(t));
    xe = zeros(size(t));
    ye = zeros(size(t));
    for i = 1:n_seg_nodes(j)-1
        % Cubic Spline
        x = (t > ti & t <= (ti + dt(j))).*(ax(i) + bx(i)*(t-ti) + cx(i)*(t-ti).^2 + dx(i)*(t-ti).^3)+x;
        y = (t > ti & t <= (ti + dt(j))).*(ay(i) + by(i)*(t-ti) + cy(i)*(t-ti).^2 + dy(i)*(t-ti).^3)+y;
        
        % Linear element edge
        xe = (t > ti & t <= (ti + dt(j))).*((t-ti)*ax(i+1)/dt(j) - (t-ti-dt(j))*ax(i)/dt(j)) + xe;
        ye = (t > ti & t <= (ti + dt(j))).*((t-ti)*ay(i+1)/dt(j) - (t-ti-dt(j))*ay(i)/dt(j)) + ye;
        
        ti = ti + dt(j);
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
    plot(x,y,'-b')
    hold all
    plot(ax,ay,'o','MarkerSize',12,'Color',[0.4660    0.6740    0.1880])
    plot(xe,ye,'k')
    
    

    hold off
end

fclose(file);




hold on

file = fopen('work/eval_nodes.out');

n = textscan(file,'%d',1);
n_seg = n{1};

for j = 1:n_seg
%     figure(j)
%     hold on
    
    n = textscan(file, '%d',1); 
    
    xy = textscan(file, '%f%f',n{1});
    
    plot(xy{:,1},xy{:,2},'x','MarkerSize',12,'Color',[0.4660    0.6740    0.1880])
end
fclose(file);











% hold on
% 
% fid = fopen('work/nodes.out') ;
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

