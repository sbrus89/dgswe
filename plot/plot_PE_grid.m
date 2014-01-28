clear all
close all
clc
restoredefaultpath

direc = '/home/sbrus/Dropbox/dgprep/work/';

lcolor = lines;

elem = 'on';
node = 'off';

for i = 1:8
    
    mesh_name = [direc,'PE000',num2str(i-1),'/fort.14'];
    
    % plotmesh
    [EToV,VX,B] = readfort14(mesh_name) ;
    p = 1 ;
    
    DEToV = zeros(length(EToV(:,1)),4) ;
    DEToV(:,1) = 3 ;
    DEToV(:,2:4) = EToV ;

    figure(1);
    hold on
    drawNGonMesh4( VX, DEToV, 'ElNum', elem, 'NodeNum', node, lcolor(i,:) );
    axis image
    ylimits = get(gca,'ylim');
    xlimits = get(gca,'xlim');

    figure(i+1);
    hold on
    drawNGonMesh4( VX, DEToV, 'ElNum', elem, 'NodeNum', node, lcolor(i,:) );
    axis image
    
end

for i = 1:8
    figure(i+1);
    ylim(ylimits);
    xlim(xlimits);
end