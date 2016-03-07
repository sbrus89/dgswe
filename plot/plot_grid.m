clear all
clc
restoredefaultpath


elcolor = 'b';
ndcolor = 'r';
lcolor = 'g';

elem = 'off';
node = 'off';

p = 1 ;

grd_direc = '/home/sbrus/Codes/dgswe/grids/';
grd_name = 'converge1_dble.grd';
% grd_name = 'inlet1.grd';

% grd_direc = '/home/sbrus/data-drive/galveston/dgswe/tri/';
% grd_name = 'galveston_tri.grd';

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,nc] = size(EToV);
[nn,~] = size(VX);


DEToV = zeros(length(EToV(:,1)),nc+1) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:nc+1) = EToV ;
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
axis image

% x1 = 2000.00000000000;
% xm = 1949.72923397938;
% x2 = 1930.52751000000;
% y1 = 400.000000000000;
% ym = 409.897987897487;
% y2 = 413.678682854092;
% 
% 
% plot([x1 x2], [y1 y2],'g')
% plot(xm,ym,'om')
% 
% t = linspace(-10,10,100);
% 
% x = xm - (ym-y1)/(xm-x1)*t;
% y = ym + t;
% 
% plot(x,y,'r')
