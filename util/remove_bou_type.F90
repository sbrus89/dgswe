      PROGRAM remove_bou_type

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name
      USE grid_file_mod
      
      IMPLICIT NONE
      
      CHARACTER(200) :: file_in 
      CHARACTER(200) :: file_out   
      
      INTEGER :: i,j
      INTEGER :: bou
      INTEGER :: bou_type
      INTEGER :: nbnds
      INTEGER, PARAMETER :: cbou_type = 30

!       file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/SMS_output/coarse.grd"
!       file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/unmodified_no_type30/coarse.grd" 
      file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/SMS_output/galveston_SL18_cart.grd"
      file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v29_cart/unmodified_no_type30/galveston_SL18_cart.grd"  
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      
      
      CALL read_header(0,file_in,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      CALL print_grid_info(file_in,grid_name,ne,nn)       
      
      
      CALL write_header(file_out,grid_name,ne,nn)
      CALL write_coords(nn,xy,depth)
      CALL write_connectivity(ne,ect,el_type,nverts)
      CALL write_open_boundaries(nope,neta,obseg,obnds)      
      
      bou = nbou
      nbou = 0
      nvel = 0        
      DO i = 1,bou
        
        bou_type = fbseg(2,i)
        nbnds = fbseg(1,i)          
          
        IF (bou_type == cbou_type) THEN
          CYCLE
        ENDIF
          
        nbou = nbou + 1         
          
        DO j = 1,nbnds
          nvel = nvel + 1
          fbnds(j,nbou) = fbnds(j,i)
        ENDDO
        fbseg(1,nbou) = nbnds
        fbseg(2,nbou) = bou_type
                    
      ENDDO
       
      CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)       
      CALL copy_footer(file_in,file_out)      
      
      
      END PROGRAM remove_bou_type