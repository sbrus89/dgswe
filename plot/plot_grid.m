clear all
clc
restoredefaultpath


elcolor = 'b';
ndcolor = 'r';
lcolor = 'r';

elem = 'off';
node = 'off';

p = 1 ;

grd_direc = '~/data-drive/converge_quad/';
grd_name = 'converge_quad4.grd';

[EToV,VX,HB,nelnds,~,~,~] = readfort14([grd_direc,grd_name]);
[ne,nc] = size(EToV);
[nn,~] = size(VX);


DEToV = zeros(length(EToV(:,1)),nc+1) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:nc+1) = EToV ;
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
axis image

