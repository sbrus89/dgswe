      PROGRAM convert2cartesian

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,r_earth,deg2rad
      USE grid_file_mod
      USE spherical_mod
      
      IMPLICIT NONE
      
      CHARACTER(100) :: grid_file
      INTEGER :: coord_sys
      REAL(rp) :: slam0,sphi0
      
      grid_file = "/home/sbrus/data-drive/galveston/SL18+TX33_andika/galveston_SL18+TX33.grd"
      coord_sys = 2
      slam0 = -94.8d0*deg2rad
      sphi0 = 29.1d0*deg2rad
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4

      CALL read_header(grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth)
      
      CALL cpp_transformation(coord_sys,r_earth,slam0,sphi0,nn,xy)

      CALL read_connectivity(ne,ect,el_type)                
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)      
      
      
      
      
      
      grid_file = "./galveston_SL18+TX33_cart.grd"
      
      CALL write_header(grid_file,grid_name,ne,nn)  
      
      CALL write_coords(nn,xy,depth)

      CALL write_connectivity(ne,ect,el_type,nverts)                
      
      CALL write_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)         
      
  

      END PROGRAM