clc

path(path,'/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/')
path(path,'/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/distmesh')
path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

% finname = '/home/sbrus/data-drive/dgswe_inlet_bath/inlet3.grd' ;
% foutname = '/home/sbrus/data-drive/dgswe_inlet_bath/inlet4.grd' ;

% finname = '/home/sbrus/Codes/dgswe/grids/converge5_dble.grd' ;
% foutname = '/home/sbrus/Codes/dgswe/grids/converge6_dble.grd' ;

finname = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/unmodified_refinements/galveston_tri.grd' 
foutname = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/unmodified_refinements/galveston_tri_x4.grd' 

% finname = '/home/sbrus/data-drive/galveston_spline_oob/grids/spline_only_refinements/galveston_tri_x4.grd' ;
% foutname = '/home/sbrus/data-drive/galveston_spline_oob/grids/unmodified_refinements/galveston_tri_x16.grd' ;

% finname = '/home/sbrus/data-drive/galveston_spline_oob/grids/spline_only_refinements/galveston_tri_x16.grd' ;
% foutname = '/home/sbrus/data-drive/galveston_spline_oob/grids/unmodified_refinements/galveston_tri_x64.grd' ;



refinefort14_alpha1( finname, foutname ) ;
% writefort14_incprec( finname, foutname ) ;


rmpath '/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/'
rmpath '/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/distmesh'
rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'
