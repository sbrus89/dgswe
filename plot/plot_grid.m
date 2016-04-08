clear all
clc
restoredefaultpath


elcolor = 'b';
ndcolor = 'r';
lcolor = 'k';

elem = 'on';
node = 'on';

p = 1 ;

% grd_direc = '/home/sbrus/Codes/dgswe/grids/';
% grd_name = 'converge1_dble.grd';
% % grd_name = 'inlet1.grd';

grd_direc = '/home/sbrus/data-drive/galveston/dgswe/tri/';
grd_name = 'galveston_tri.grd';

% grd_direc = '/home/sbrus/data-drive/galveston/adcirc/refinedx64/ESL0/';
% grd_name = 'fort.14';

[EToV,VX,HB,nelnds,opedat,boudat,title] = readfort14([grd_direc,grd_name]);
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






% figure
% hold on
% for bou = 1:boudat.nbou
%    for nd = 1:boudat.nvell(bou)-1
%        nd1 = boudat.nbvv(nd,bou);
%        nd2 = boudat.nbvv(nd+1,bou);
%        plot([VX(nd1,1),VX(nd2,1)], [VX(nd1,2),VX(nd2,2)])
%    end
%     
% end
% axis equal