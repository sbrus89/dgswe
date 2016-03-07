clear all
close all
clc

makemovie = 0;
nsnap = 100;
%direc = '../curvilinear/';
direc = '../ptr_test2/';
% direc = '../quad_loop/';
%direc = '../LDG/';
%direc = '../openmp/';

% direc = '..';

FramesFolder = strcat(direc,'output/velsnap') ;
if ( exist(FramesFolder,'dir') == 0 ) 
    mkdir(FramesFolder) ; 
end

% f1 = figure('Position',[1360 667 560 420]);
% f2 = figure('Position',[798 667 560 420]);

h = figure('Position',[100 100 1500 1000]) ;

% Open water column height, x momentum, and y momentum modal solution files
fid_H = fopen([direc,'output/solution_H.d']);
fid_Qx = fopen([direc,'output/solution_Qx.d']);
fid_Qy = fopen([direc,'output/solution_Qy.d']);

% Read in grid name
name = fscanf(fid_H,' %s ', 1);
name = fscanf(fid_Qx,' %s ', 1);
name = fscanf(fid_Qy,' %s ', 1);

% Read in grid file
[EToV,VX,HB,~,~,~] = readfort14([direc,name]);

% Open L2 projection information file
fid = fopen([direc,'output/projection.d']);

N = fscanf(fid, '%g', 1); % number of degrees of freedom
mL2 = fscanf(fid,' %g ',[N 3])'; % L2 projection matrix mL2(i,j) = integral of (phi_DG_basis(j)*phi_linear_nodal(i))
mm = fscanf(fid,' %g ',[3 3])'; % linear nodal mass matrix
fclose(fid);

[ne,~] = size(EToV); % number of elements
[nn,~] = size(VX); % number of nodes

xmin = min(VX(:,1));
xmax = max(VX(:,1));

ymin = min(VX(:,2));
ymax = max(VX(:,2));

nx = 300;
ny = 300;

[X,Y] = meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

area = zeros(ne,1);
for el = 1:ne
    x1 = VX(EToV(el,1),1);
    x2 = VX(EToV(el,2),1);
    x3 = VX(EToV(el,3),1);
    
    y1 = VX(EToV(el,1),2);
    y2 = VX(EToV(el,2),2);
    y3 = VX(EToV(el,3),2);
    
    area(el) = .5*((x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1));
end

nepn = zeros(1,nn);
for el = 1:ne
    n1 = EToV(el,1);
    n2 = EToV(el,2);
    n3 = EToV(el,3);
    
    nepn(n1) = nepn(n1) + 1;   
    nepn(n2) = nepn(n2) + 1;    
    nepn(n3) = nepn(n3) + 1;
end
mnepn = max(nepn);

nepn = zeros(1,nn);
epn = zeros(mnepn,nn);
for el = 1:ne
    n1 = EToV(el,1);
    n2 = EToV(el,2);
    n3 = EToV(el,3);
    
    nepn(n1) = nepn(n1) + 1;
    epn(nepn(n1),n1) = el;
    
    nepn(n2) = nepn(n2) + 1;
    epn(nepn(n2),n2) = el;
    
    nepn(n3) = nepn(n3) + 1;
    epn(nepn(n3),n3) = el;
    
end

kdtreeobj = KDTreeSearcher(VX);

R = zeros(size(X));
S = zeros(size(Y));
elin = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        xpt = X(i,j);
        ypt = Y(i,j);

        n = knnsearch(kdtreeobj,[xpt ypt]);
        
        for el = 1:nepn(n)
            indomain = 0;
            
            x1 = VX(EToV(epn(el,n),1),1);
            x2 = VX(EToV(epn(el,n),2),1);
            x3 = VX(EToV(epn(el,n),3),1);
            
            y1 = VX(EToV(epn(el,n),1),2);
            y2 = VX(EToV(epn(el,n),2),2);
            y3 = VX(EToV(epn(el,n),3),2);
            
            xp = [x1 x2 x3];
            yp = [y1 y2 y3];
            
            [in,on] = inpolygon(xpt,ypt,xp,yp);
            
            if in == 1 || on == 1
               indomain = 1;
               break
            end
        end
        
        if indomain == 1
            r = ((y3-y1)*(xpt-.5*(x2+x3))+(x1-x3)*(ypt-.5*(y2+y3)))/area(epn(el,n));
            s = ((y1-y2)*(xpt-.5*(x2+x3))+(x2-x1)*(ypt-.5*(y2+y3)))/area(epn(el,n));
            
            R(i,j) = r;
            S(i,j) = s;
            elin(i,j) = epn(el,n);
        else
            elin(i,j) = ne + 1;
        end
    end
end

Hsol = zeros(nx,ny,nsnap);
u = zeros(nx,ny,nsnap);
v = zeros(nx,ny,nsnap);
t = zeros(nsnap);

snap = 0;
while ~feof(fid_H) && snap < nsnap

snap = snap + 1;
    
th = fscanf(fid_H,' %g ', 1); % read in time
H = fscanf(fid_H,' %g ', [ne N])'; % read in H solution at time t

tqx = fscanf(fid_Qx,' %g ', 1); % read in time
Qx = fscanf(fid_Qx,' %g ', [ne N])'; % read in Qx solution at time t

t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
Qy = fscanf(fid_Qy,' %g ', [ne N])'; % read in Qy solution at time t

for i = 1:nx
    for j = 1:ny
        
        if elin(i,j) <= ne
            phi = basis_eval(R(i,j),S(i,j));
            
            Hsol(i,j,snap) = phi(1,1:N)*H(:,elin(i,j));
            u(i,j,snap) = (phi(1,1:N)*Qx(:,elin(i,j)))/Hsol(i,j);
            v(i,j,snap) = (phi(1,1:N)*Qy(:,elin(i,j)))/Hsol(i,j);
        else
            Hsol(i,j,snap) = NaN;
            u(i,j,snap) = NaN;
            v(i,j,snap) = NaN;
        end
    end
end


end

vel = sqrt(u.^2 + v.^2);

% projsol = zeros(nn,4,nsnap);
% projsol(:,1,:) = zeta;
% projsol(:,2,:) = u;
% projsol(:,3,:) = v;
% projsol(:,4,:) = vel;

for tsnap = 1:snap
    
    disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    [day,hr,min,sec] = s2dhms(t(tsnap));
    

    contourf(X,Y,vel(:,:,tsnap),50,'LineStyle','none') ;
    colorbar;
%     hold on
%     quiver( X, Y,u(:,:,tsnap), v(:,:,tsnap), '-k' ) ;
%     grid on
%     hold off
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(min),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    axis image
    %zlim([-.5 1])
    %caxis([0 1])

    if makemovie == 1
       if tsnap == 1
           set(gca,'nextplot','replacechildren','visible','off')
            f = getframe(2);
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,20) = 0;
       else
            f = getframe(2);
            im(:,:,1,tsnap) = rgb2ind(f.cdata,map,'nodither');
       end
    end

    pause(.01)

    print('-r300','-dpng',sprintf('%s/surf%04d',FramesFolder,tsnap)) ;
    
end

if makemovie == 1
    imwrite(im,map,'Inlet.gif','DelayTime',0,'LoopCount',inf)
end
