clear all
clc
restoredefaultpath

path(path,'/home/sbrus/Dropbox/Code/2Dswe_tri/adcirc/inlet/mesh1/p1')
% direc = '/home/sbrus/Dropbox/dgswe/output/PE0001';
% path(path,direc)


lcolor = 'r';

elem = 'off';
node = 'on';

% plotmesh
[EToV,VX,B] = readfort14() ;
p = 1 ;

% draw mesh
% [G,ElIdDG,xDG,vKN,vKI] = MakeNodalDGNodes( p, EToV, VX ) ;

DEToV = zeros(length(EToV(:,1)),4) ;
DEToV(:,1) = 3 ;
DEToV(:,2:4) = EToV ;
drawNGonMesh4( VX, DEToV, 'ElNum', elem, 'NodeNum', node, lcolor )
axis image

rmpath '/home/sbrus/Dropbox/Code/2Dswe_tri/adcirc/inlet/mesh1/p1'
% rmpath direc


