      SUBROUTINE read_grid()
      
      USE globals, ONLY: rp,ne,nn,ect,vct,xy,depth,nelnds,elxy,elhb,bndxy, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,mnelnds,mnnds, &
                         r_earth,deg2rad
                             
      USE messenger2, ONLY: myrank
      USE read_dginp, ONLY: grid_file,bathy_file,curve_file, &
                            cb_file_exists,hb_file_exists, &
                            ctp,hbp,h0,coord_sys,slam0,sphi0
                            
      USE grid_file_mod                            

      IMPLICIT NONE




      CALL read_header(grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth,h0,coord_sys,r_earth,slam0,sphi0)

      CALL read_connectivity(ne,ect,el_type,nelnds,mnelnds)
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      CALL read_bathy_file(myrank,bathy_file,hbp,ne,elhb,depth,ect,nelnds)
                  
      CALL read_curve_file(myrank,curve_file,ctp,nbou,xy,bndxy)      
      
      IF (myrank == 0) THEN
        CALL print_grid_info(grid_file,grid_name,ne,nn)
      ENDIF
      
      
      DO i = 1,ne        
        DO j = 1,nelnds(el)
          elxy(j,el,1) = xy(1,ect(j,el))
          elxy(j,el,2) = xy(2,ect(j,el))
        ENDDO              
      ENDDO       
      

      RETURN
      END SUBROUTINE read_grid