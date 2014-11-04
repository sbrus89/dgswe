clear all
close all
clc
restoredefaultpath

direc = '/home/sbrus/Codes/dgswe/work/';

lcolor = lines;

elem = 'on';
elcolor = 'r';
node = 'off';
ndcolor = 'k';

npe = 0;
found = 1;
while found == 1    
    if exist([direc,'PE000',num2str(npe)],'dir')
        npe = npe + 1;
    else
        found = 0;
    end
end

for i = 1:npe
    
    mesh_name = [direc,'PE000',num2str(i-1),'/fort.14'];
    
    % plotmesh
    [EToV,VX,B,nelnds] = readfort14(mesh_name) ;
    p = 1 ;
    
    DEToV = zeros(length(EToV),10);
    DEToV(:,1) = nelnds;
    DEToV(:,2:end) = EToV ;

    figure(1);
    hold on
    drawNGonMesh4( VX, DEToV, lcolor(i,:), 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )  ;  
    axis image
    ylimits = get(gca,'ylim');
    xlimits = get(gca,'xlim');

    figure(i+1);
    hold on
    drawNGonMesh4( VX, DEToV, lcolor(i,:), 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )  ; 
    axis image
    
end

for i = 1:npe
    figure(i+1);
    ylim(ylimits);
    xlim(xlimits);
end