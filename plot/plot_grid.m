clear all
clc
restoredefaultpath


elcolor = 'b';
ndcolor = 'r';
lcolor = 'k';

elem = 'on';
node = 'on';

p = 1 ;

grd_direc = '/home/sbrus/data-drive/dgswe_inlet_bath/';
% grd_direc = '/home/sbrus/Codes/dgswe/grids/';
grd_name = 'inlet1.grd';

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,nc] = size(EToV);
[nn,~] = size(VX);


DEToV = zeros(length(EToV(:,1)),nc+1) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:nc+1) = EToV ;
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
axis image

