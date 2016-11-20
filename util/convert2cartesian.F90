      PROGRAM convert2cartesian

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,r_earth,deg2rad
      USE grid_file_mod
      USE spherical_mod
      
      IMPLICIT NONE
      
      INTEGER :: i
      CHARACTER(100) :: grid_file
      INTEGER :: coord_sys
      REAL(rp) :: slam0,sphi0
      REAL(rp) :: a,b,r

      
      grid_file = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v19/galveston_SL18.grd"

      coord_sys = 2
      slam0 = -94.90d0*deg2rad
      sphi0 = 29.195d0*deg2rad
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4

      CALL read_header(0,grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth)
      
      CALL cpp_transformation(coord_sys,r_earth,slam0,sphi0,nn,xy)

!       a = 1.708410418613708d5
!       b = 1.651273587295530d7
!       r = 1.106591847107369d5
!       DO i = 1,nn
!         xy(1,i) = a*xy(1,i) + b
!         xy(2,i) = r*xy(2,i)
!       ENDDO

      CALL read_connectivity(ne,ect,el_type)                
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)      
      
      
      
      
      
      grid_file = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v19_cart/galveston_SL18_cart.grd"
      
!       CALL cpp_inverse(coord_sys,r_earth,slam0,sphi0,nn,xy)      
      
      CALL write_header(grid_file,grid_name,ne,nn)  
      
      CALL write_coords(nn,xy,depth)

      CALL write_connectivity(ne,ect,el_type,nverts)                
      
      CALL write_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)         
      
  

      END PROGRAM