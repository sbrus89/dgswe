clear all
close all
clc
restoredefaultpath

% direc = '/home/sbrus/Codes/dgswe/work/';
% direc = '/home/sbrus/data-drive/galveston_spline/tri/p3/ctp3/hbp3/test/';
% direc = '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x4/p3/ctp3/hbp3/';
% direc = '/home/sbrus/data-drive/galveston_spline/tri_x64/p3/ctp3/hbp3/';
% direc = '/home/sbrus/data-drive/EC2001/dgswe/';
direc = '/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/rimls_spline_modified/';

lcolor = lines;

elem = 'off';
elcolor = 'k';
node = 'off';
ndcolor = 'k';

npe = 0;
found = 1;
while found == 1    
    pe = sprintf('%04d',npe);
    if exist([direc,'PE',pe],'dir')
        npe = npe + 1;
    else
        found = 0;
    end
end

for i = 1:npe
    
    mesh_name = [direc,'PE',sprintf('%04d',i-1),'/fort.14'];
    
    cl = mod(i,length(lcolor))+1;
    
    % plotmesh
    [EToV,VX,B,nelnds,~,~,~] = readfort14(mesh_name) ;
    p = 1 ;
    
    [nel,~] = size(EToV);
    DEToV = zeros(nel,10);
    DEToV(:,1) = nelnds;
    DEToV(:,2:end) = EToV ;

    figure(1);
    hold on
    drawNGonMesh4( VX, DEToV, lcolor(cl,:), 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )  ;  
    axis equal
    
    xmax = max(VX(:,1));
    xmin = min(VX(:,1));
    ymax = max(VX(:,2));
    ymin = min(VX(:,2));
    
    text(.5*(xmin+xmax),.5*(ymin+ymax),num2str(i-1));
    
    
    ylimits = get(gca,'ylim');
    xlimits = get(gca,'xlim');

%     figure(i+1);
%     hold on
%     drawNGonMesh4( VX, DEToV, lcolor(i,:), 'ElNum', elem, elcolor, 'NodeNum', node, ndcolor )  ; 
%     axis image
    
end



% for i = 1:npe
%     figure(i+1);
%     ylim(ylimits);
%     xlim(xlimits);
% end