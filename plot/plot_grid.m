clear all
clc
restoredefaultpath


elcolor = 'r';
ndcolor = 'g';
lcolor = 'k';

elem = 'on';
node = 'on';

p = 1 ;

% grd_direc = '/home/sbrus/Codes/dgswe/grids/';
% grd_name = 'converge5_dble.grd';
% grd_name = 'inlet1.grd';

% grd_direc = '/home/sbrus/data-drive/galveston_spline/grids/';
% grd_name = 'galveston_tri_x64.grd';

% grd_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/';
% grd_name = 'fort.14_refine';

% grd_direc = '/home/sbrus/data-drive/galveston_spline_oob/grids/unmodified_refinements/';
% grd_name = 'galveston_tri_x4.grd';

% grd_direc = '/home/sbrus/data-drive/galveston_spline_oob/grids/spline_only_refinements/';
% grd_name = 'galveston_tri_x64.grd';

% grd_direc = '/home/sbrus/data-drive/galveston_SL18/';
% grd_name = 'galveston_SL18.grd';

%  grd_direc ='/home/sbrus/data-drive/galveston/dgswe/quad2/';
%  grd_name = 'galveston_quad2.grd';
 
% grd_direc = '/home/sbrus/data-drive/galveston_spline/grids/rimls/';
% grd_name = 'fort.14_rimls';

% grd_direc = '/home/sbrus/Codes/dgswe/spline/work/';
% grd_name = 'nodes.out';

grd_direc = '/home/sbrus/Codes/dgswe/work/PE0182/';
grd_name = 'fort.14';

[EToV,VX,HB,nelnds,opedat,boudat,title] = readfort14([grd_direc,grd_name]);
[ne,nc] = size(EToV);
[nn,~] = size(VX);


DEToV = zeros(length(EToV(:,1)),nc+1) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:nc+1) = EToV ;
% drawNGonMesh4_bou( VX, DEToV,boudat.nbvv, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
axis equal

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