clear 
clc
close all

path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')

[EToV,VX,B] = readfort14('/home/sbrus/Codes/dgswe/grids/inlet1.grd') ;

fid = fopen('maxdiff_el.d');
ne = fscanf(fid,'%g',1) ;
% diff = zeros(ne,3) ;
diff = fscanf(fid,'%g %g %g\n', [3 ne])' ;

figure(1)
hold on
for el = 1:ne    
  fill(VX(EToV(el,:),1),VX(EToV(el,:),2),diff(el,1))   
end
colorbar
axis image

figure(2)
hold on
for el = 1:ne    
  fill(VX(EToV(el,:),1),VX(EToV(el,:),2),diff(el,2))   
end
colorbar
caxis([0 10^-12])
axis image

figure(3)
hold on
for el = 1:ne    
  fill(VX(EToV(el,:),1),VX(EToV(el,:),2),diff(el,3))   
end
colorbar
axis image

rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'