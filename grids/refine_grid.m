path(path,'/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/')
path(path,'/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/distmesh')
path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

finname = '/home/sbrus/data-drive/dgswe_inlet_bath/inlet3.grd' ;
foutname = '/home/sbrus/data-drive/dgswe_inlet_bath/inlet4.grd' ;

refinefort14_alpha1( finname, foutname ) ;
% writefort14_incprec( finname, foutname ) ;


rmpath '/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/'
rmpath '/home/sbrus/Codes/SomeDGMaterials/MeshUTIL/distmesh'
rmpath '/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG'
rmpath '/home/sbrus/Codes/SomeDGMaterials/NodalDG'