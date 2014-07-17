clear all
clc
restoredefaultpath


elcolor = 'b';
ndcolor = 'r';
lcolor = 'r';

elem = 'off';
node = 'off';

% plotmesh
% [EToV,VX,B] = read_fort14() ;
p = 1 ;

% filename = '~/tri2quad/grids/converge2_quad.grd';
filename = '~/converge_quad/adjust_grid/work/converge2_quad.grd';

fid = fopen(filename);

agrid = fgetl(fid) ;
disp(agrid) ;
title = agrid ; 

N = fscanf(fid,'%g %g',2) ;

Val = zeros(N(2),4) ;
% Nov 15, 2012, improve reading efficiency
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;

iv = sort(Val(:,1)) ;
Val = Val(iv,:) ; 

idx = zeros(N(1),6) ; 
% Nov 15, 2012, improve reading efficient
idx = fscanf(fid,'%d %d %d %d %d %d\n', [6 N(1)])' ;
iv = sort(idx(:,1)) ;
idx = idx(iv,:) ; 

VX = zeros(N(2),2) ;
B = zeros(N(2),1) ;
EToV = zeros(N(1),4) ;

% Arrange it to a Nodal DG input
VX = Val(:,2:3) ;
B  = Val(:,4) ;
EToV = idx(:,3:6) ; 
nelnds = idx(:,2);


DEToV = zeros(length(EToV(:,1)),5) ;
DEToV(:,1) = nelnds ;
DEToV(:,2:5) = EToV ;
drawNGonMesh4( VX, DEToV, lcolor, 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )
axis image

