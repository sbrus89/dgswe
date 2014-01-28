close all
clear all
clc

t = 0:.00001:1;

file = fopen('../curvilinear/spline.out');
file2 = fopen('../curvilinear/norm.out');
file3 = fopen('../curvilinear/norm2.out');

figure

n = textscan(file,'%d',1);
n_seg = n{1};

n = textscan(file2,'%d%d',1);
nqpte = n{2};

n = textscan(file3,'%d%d',1);
nqpte = n{2};

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
   
   norms = textscan(file2, '%f%f%f%f%f%f%f',(n_seg_nodes(j)-1)*nqpte);
   tp = norms{:,1};
   nx = norms{:,2};
   ny = norms{:,3};
   nxedge = norms{:,4};
   nyedge = norms{:,5};
   nlx = norms{:,6};
   nly = norms{:,7};
    
   norms = textscan(file3, '%f%f%f%f',(n_seg_nodes(j)-1)*nqpte);
   xp = norms{:,1};
   yp = norms{:,2}; 
   nx2 = norms{:,3};
   ny2 = norms{:,4};
   
   %%%%%%%%%% Plot boundary polynomials %%%%%%%%%%
    ti = 0;
    x = zeros(size(t));
    y = zeros(size(t));
    xe = zeros(size(t));
    ye = zeros(size(t));
    xl = zeros(size(t));
    yl = zeros(size(t));
    for i = 1:n_seg_nodes(j)-1
        % Cubic Spline
        x = (t > ti & t <= (ti + dt(j))).*(ax(i) + bx(i)*(t-ti) + cx(i)*(t-ti).^2 + dx(i)*(t-ti).^3)+x;
        y = (t > ti & t <= (ti + dt(j))).*(ay(i) + by(i)*(t-ti) + cy(i)*(t-ti).^2 + dy(i)*(t-ti).^3)+y;
        
        % Linear element edge
        xe = (t > ti & t <= (ti + dt(j))).*((t-ti)*ax(i+1)/dt(j) - (t-ti-dt(j))*ax(i)/dt(j)) + xe;
        ye = (t > ti & t <= (ti + dt(j))).*((t-ti)*ay(i+1)/dt(j) - (t-ti-dt(j))*ay(i)/dt(j)) + ye;
        
%         % Lagange Cubic
%         dx1 = ax(i+1) - ax(i);
%         dx2 = ax(i+2) - 2*ax(i+1) + ax(i);
%         dx3 = ax(i+3) - 3*ax(i+2) + 3*ax(i+1) - ax(i);
%         
%         dy1 = ay(i+1) - ay(i);
%         dy2 = ay(i+2) - 2*ay(i+1) + ay(i);
%         dy3 = ay(i+3) - 3*ay(i+2) + 3*ay(i+1) - ay(i);
%         
%         xl = (t > ti & t <= (ti + dt(j))).*(ax(i) + (t-ti)*dx1/dt(j) + (t-ti).*(t-(ti+dt(j)))*dx2/(2*dt(j)^2)+ (t-ti).*(t-(ti+dt(j))).*(t-(ti+2*dt(j)))*dx3/(6*dt(j)^3)) + xl;
%         yl = (t > ti & t <= (ti + dt(j))).*(ay(i) + (t-ti)*dy1/dt(j) + (t-ti).*(t-(ti+dt(j)))*dy2/(2*dt(j)^2)+ (t-ti).*(t-(ti+dt(j))).*(t-(ti+2*dt(j)))*dy3/(6*dt(j)^3)) + yl;
        
        ti = ti + dt(j);
    end
    hold on
    plot(x(2:end-1),y(2:end-1))
%    plot(xl(2:end-1),yl(2:end-1))
    hold all
    plot(ax,ay,'rx')
    plot(xe,ye,'k')
    
    
    %%%%%%%% Plot normals %%%%%%%%
    hold on
    ti = 0;
    m = 1;
    for i = 1:n_seg_nodes(j)-1
        for k = 1:nqpte
            
            % Cubic spline normals
            xpt(1) = ax(i) + bx(i)*(tp(m)-ti) + cx(i)*(tp(m)-ti).^2 + dx(i)*(tp(m)-ti).^3 ;
            ypt(1) = ay(i) + by(i)*(tp(m)-ti) + cy(i)*(tp(m)-ti).^2 + dy(i)*(tp(m)-ti).^3 ;
        
            xpt(2) = xpt(1) + 50*nx(m);
            ypt(2) = ypt(1) + 50*ny(m);
            
            xpt2(1) = xp(m);
            ypt2(1) = yp(m);
            
            xpt2(2) = xpt2(1) + 50*nx2(m);
            ypt2(2) = ypt2(1) + 50*ny2(m);
            
            % Linear edge normals
            xpedge(1) = ((tp(m)-ti)*ax(i+1)/dt(j) - (tp(m)-ti-dt(j))*ax(i)/dt(j));
            ypedge(1) = ((tp(m)-ti)*ay(i+1)/dt(j) - (tp(m)-ti-dt(j))*ay(i)/dt(j));
            
            xpedge(2) = xpedge(1) + 50*nxedge(m);
            ypedge(2) = ypedge(1) + 50*nyedge(m);
            
%             % Cubic Lagrange normals
%             dx1 = ax(i+1) - ax(i);
%             dx2 = ax(i+2) - 2*ax(i+1) + ax(i);
%             dx3 = ax(i+3) - 3*ax(i+2) + 3*ax(i+1) - ax(i);
%         
%             dy1 = ay(i+1) - ay(i);
%             dy2 = ay(i+2) - 2*ay(i+1) + ay(i);
%             dy3 = ay(i+3) - 3*ay(i+2) + 3*ay(i+1) - ay(i);
%             
%             xlp(1) = (ax(i) + (tp(m)-ti)*dx1/dt(j) + (tp(m)-ti).*(tp(m)-(ti+dt(j)))*dx2/(2*dt(j)^2) + (tp(m)-ti).*(tp(m)-(ti+dt(j))).*(tp(m)-(ti+2*dt(j)))*dx3/(6*dt(j)^3));
%             ylp(1) = (ay(i) + (tp(m)-ti)*dy1/dt(j) + (tp(m)-ti).*(tp(m)-(ti+dt(j)))*dy2/(2*dt(j)^2) + (tp(m)-ti).*(tp(m)-(ti+dt(j))).*(tp(m)-(ti+2*dt(j)))*dy3/(6*dt(j)^3));
%             
%             xlp(2) = xlp(1) + 50*nlx(m);
%             ylp(2) = ylp(1) + 50*nly(m);
        
            plot(xpt,ypt,'g')
            plot(xpt2,ypt2,'m')
            plot(xpedge,ypedge,'k')
      %      plot(xlp,ylp,'m')
            
            m = m+1;
        end
        ti = ti + dt(j);
    end


    hold off
end

fclose(file);
fclose(file2);
fclose(file3);


axis equal
%axis([23000 27000 10000 12500])

hold on
load ../curvilinear/coords.out
plot(coords(:,1),coords(:,2),'x')
plotfort14


