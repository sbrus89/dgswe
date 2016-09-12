      PROGRAM convert_open_boundaries

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name
      USE grid_file_mod

      IMPLICIT NONE
      
      INTEGER :: i,j
      INTEGER :: nbou_mod,nvel_mod
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbseg_mod
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbnds_mod
      CHARACTER(100) :: grid_file_in,grid_file_out
      
      grid_file_in = "/home/sbrus/data-drive/galveston_spline_flux/grids/spline_mls_modified/galveston_tri_x16.grd"
      grid_file_out = "/home/sbrus/data-drive/galveston_spline_flux/grids/galveston_tri_x16.grd"
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4

      CALL read_header(0,grid_file_in,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      

      
      ALLOCATE(fbseg_mod(2,nbou+1))
      ALLOCATE(fbnds_mod(nvel+neta,nbou+1))
      
      DO i = 1,nbou
        fbseg_mod(1,i) = fbseg(1,i)
        fbseg_mod(2,i) = fbseg(2,i)
        DO j = 1,fbseg(1,i)
          fbnds_mod(j,i) = fbnds(j,i)          
        ENDDO
      ENDDO
      
      nbou_mod = nbou + 1
      nvel_mod = nvel + neta   
      
      nope = 0
      neta = 0
      
      fbseg_mod(1,nbou_mod) = obseg(1) 
      fbseg_mod(2,nbou_mod) = 12
      DO j = 1,obseg(1)
        fbnds_mod(j,nbou_mod) = obnds(j,1)
      ENDDO
      
      CALL write_header(grid_file_out,grid_name,ne,nn)
      CALL write_coords(nn,xy,depth)
      CALL write_connectivity(ne,ect,el_type,nverts)
      CALL write_open_boundaries(nope,neta,obseg,obnds)
      CALL write_flow_boundaries(nbou_mod,nvel_mod,fbseg_mod,fbnds_mod)

      END PROGRAM convert_open_boundaries