clear all
clc
restoredefaultpath


elcolor = 'r';
ndcolor = 'r';
lcolor = 'r';

elem = 'on';
node = 'off';

p = 1 ;

% grd_direc = '/home/sbrus/data-drive/dgswe_inlet_bath/';
grd_direc = '/home/sbrus/Codes/dgswe/grids/';
grd_name = 'converge2_dble.grd';

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,nc] = size(EToV);
[nn,~] = size(VX);


DEToV = zeros(length(EToV(:,1)),nc+1) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:nc+1) = EToV ;
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
% axis image

