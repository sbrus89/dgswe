      PROGRAM adjust

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name
                         
      USE grid_file_mod

      IMPLICIT NONE
      
      INTEGER :: nd,bnd
      INTEGER :: nbnds,btype,n
      REAL(rp) :: x,y,ytest
      CHARACTER(100) :: grid_file
      
      
      grid_file = "/home/sbrus/Codes/dgswe/grids/converge6_dble.grd"
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4

      CALL read_header(grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth)      

      CALL read_connectivity(ne,ect,el_type)                
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)      
      
      
      DO bnd = 1,nbou
        nbnds = fbseg(1,bnd)
        btype = fbseg(2,bnd)
        
        IF (btype == 0 .OR. btype == 10 .OR. btype == 20) THEN
          DO nd = 1,nbnds
            n = fbnds(nd,bnd)
          
            x = xy(1,n)
            ytest = xy(2,n)
          
            IF (ytest < 250d0) THEN
              y = 0d0 + 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
            ELSE IF (ytest > 250d0) THEN
              y = 500d0 - 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
            ENDIF
          
            xy(2,n) = y
          
          ENDDO
        ENDIF
      ENDDO
      
      
      CALL write_header(grid_file,grid_name,ne,nn)  
      
      CALL write_coords(nn,xy,depth)

      CALL write_connectivity(ne,ect,el_type,nverts)                
      
      CALL write_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)  

      END PROGRAM adjust
      
      