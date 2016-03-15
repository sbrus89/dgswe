      SUBROUTINE read_grid()
      
      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,elxy,elhb,bndxy, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,mnnds, &
                         r_earth,deg2rad
                             
      USE messenger2, ONLY: myrank
      USE read_dginp, ONLY: grid_file,bathy_file,curve_file, &
                            cb_file_exists,hb_file_exists, &
                            ctp,hbp,h0,coord_sys,slam0,sphi0
                            
      USE grid_file_mod                            

      IMPLICIT NONE




      CALL read_header(grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth,h0,coord_sys,r_earth,slam0,sphi0)

      CALL read_connectivity(ne,ect,el_type)
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      CALL read_bathy_file(myrank,bathy_file,hbp,ne,el_type,nverts,depth,ect,elhb)
                  
      CALL read_curve_file(myrank,curve_file,ctp,nbou,xy,bndxy)      
      
      IF (myrank == 0) THEN
        CALL print_grid_info(grid_file,grid_name,ne,nn)
      ENDIF
      
      
  
      

      RETURN
      END SUBROUTINE read_grid