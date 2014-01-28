path(path,'/home/sbrus/Dropbox/converge-curve/MeshUTIL/')
path(path,'/home/sbrus/Dropbox/converge-curve/MeshUTIL/distmesh')
path(path,'/home/sbrus/Documents/SomeDGMaterials/ModNodeADCIRCDG')
path(path,'/home/sbrus/Documents/SomeDGMaterials/NodalDG')

finname = 'inlet2.grd' ;
foutname = 'inlet3.grd' ;
refinefort14_alpha1( finname, foutname ) ;

rmpath '/home/sbrus/Dropbox/converge-curve/MeshUTIL/'
rmpath '/home/sbrus/Dropbox/converge-curve/MeshUTIL/distmesh'
rmpath '/home/sbrus/Documents/SomeDGMaterials/ModNodeADCIRCDG'
rmpath '/home/sbrus/Documents/SomeDGMaterials/NodalDG'