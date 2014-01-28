clear all
%close all
clc

makemovie = 0;
nsnap = 50;

dir = '../adcirc/converge/';

f1 = figure('Position',[1360 667 560 420]);
f2 = figure('Position',[798 667 560 420]);

% Read in grid file
[EToV,VX,HB,~,~,~] = readfort14([dir,'fort.14']);

% Open L2 projection information file
fid = fopen([dir,'projection.d']);

N = fscanf(fid, '%g', 1); % number of degrees of freedom
mL2 = fscanf(fid,' %g ',[N 3])'; % L2 projection matrix mL2(i,j) = integral of (phi_DG_basis(j)*phi_linear_nodal(i))
mm = fscanf(fid,' %g ',[3 3])'; % linear nodal mass matrix
fclose(fid);

[ne ~] = size(EToV); % number of elements
[nn ~] = size(VX); % number of nodes

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

% Assemble global linear mass matrix
M = zeros(nn);
for el = 1:ne
    for i = 1:3
        for j = 1:3
            M(EToV(el,i),EToV(el,j)) = M(EToV(el,i),EToV(el,j)) + .5*area(el)*mm(i,j);
        end
    end
end

M = sparse(M);

[ze,tt,nsnapi,Nel,ntot] = readdg63( nsnap,dir );
[Q,tt,nsnap,Nel,ntot] = readdg64( nsnapi,dir );

H = zeros(N,ne);
Qx = zeros(N,ne);
Qy = zeros(N,ne);

zeta = zeros(nn,nsnap);
u = zeros(nn,nsnap);
v = zeros(nn,nsnap);
vel = zeros(nn,nsnap);
t = zeros(nsnap);

for snap = 1:nsnap
 
    k = 1;
    for el = 1:ne
        for i = 1:N
            H(i,el) = ze(k,1,snap);
            Qx(i,el) = Q(k,1,snap);
            Qy(i,el) = Q(k,2,snap);
            k = k+1;
        end
    end

prod_H = mL2*H; % compute product of elemental degrees of freedom and projection matrix
prod_Qx = mL2*Qx;
prod_Qy = mL2*Qy;

% Assemble RHS vector
B_H = zeros(nn,1);
B_Qx = zeros(nn,1);
B_Qy = zeros(nn,1);
for el = 1:ne
    for i = 1:3
       B_H(EToV(el,i)) = B_H(EToV(el,i)) + .5*area(el)*prod_H(i,el);
       B_Qx(EToV(el,i)) = B_Qx(EToV(el,i)) + .5*area(el)*prod_Qx(i,el); 
       B_Qy(EToV(el,i)) = B_Qy(EToV(el,i)) + .5*area(el)*prod_Qy(i,el); 
    end
end

HL2 = M\B_H; % Compute L2 projected solution
QxL2 = M\B_Qx;
QyL2 = M\B_Qy;

zeta(:,snap) = HL2;
u(:,snap) = QxL2./(HL2+HB);
v(:,snap) = QyL2./(HL2+HB);
vel(:,snap) = sqrt(u(:,snap).^2 + v(:,snap).^2);
end

for tsnap = 1:snap
    
    figure(f1)
    %pdeplot( VX', [], EToV','xydata',zeta(:,tsnap),'zdata',zeta(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
    pdeplot( VX', [], EToV','xydata',zeta(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
    grid on
    [day hr min sec] = s2dhms(tt(tsnap,1));
    ttext = ['Elevation solution: t = ',num2str(tt(tsnap,1)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(min),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    axis image
    %zlim([-.5 1])
    %caxis([3.5 5])

    figure(f2)
    pdeplot( VX', [], EToV','xydata',vel(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
    hold on
    quiver( VX(:,1), VX(:,2),u(:,snap), v(:,snap), '-k' ) ;
    grid on
    hold off
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(min),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    axis image
    %zlim([-.5 1])
    %caxis([3.5 5])

    if makemovie == 1
       if tsnap == 1
           set(gca,'nextplot','replacechildren','visible','off')
            f = getframe;
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,20) = 0;
       else
            f = getframe;
            im(:,:,1,tsnap) = rgb2ind(f.cdata,map,'nodither');
       end
    end

    pause(.01)

end

if makemovie == 1
    imwrite(im,map,'DamBreak.gif','DelayTime',0,'LoopCount',inf)
end
