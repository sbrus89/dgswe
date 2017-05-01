      PROGRAM write_bou_file

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         nobfr,obtag,obfreq,obnfact,obeq,obamp,obph, &
                         nfbfr,fbtag,fbfreq,fbnfact,fbeq,fbamp,fbph, &
                         nfbsfr,fbstag,fbsbgn,fbsend,fbssig,fbsamp
      USE grid_file_mod
      USE allocation, ONLY: alloc_forcing_arrays
      
      
      IMPLICIT NONE      
      
      CHARACTER(200) :: grid_file 
      CHARACTER(200) :: bou_file   
      
      INTEGER :: i,j
      INTEGER :: nd,bfr,bou
      INTEGER :: nfrbou
      INTEGER :: segtype
      INTEGER :: nbnds
      
!       grid_file = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/galveston_SL18_cart.grd"
! !       bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/galveston_SL18_cart.bfr"         
!       bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/galveston_SL18_cart.bfr_surge"        
      
!       grid_file = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse.grd"
! !       bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse.bfr"  
!       bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse.bfr_surge"
      
      grid_file = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse_x2.grd"
!       bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse_x2.bfr"      
      bou_file  = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v27_cart/rimls_spline_modified/coarse_x2.bfr_surge" 
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      
      
      CALL read_header(0,grid_file,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      CALL print_grid_info(grid_file,grid_name,ne,nn)        
            
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      nobfr = 0
      CALL alloc_forcing_arrays(1)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
      
      nfbfr = 1
      CALL alloc_forcing_arrays(2)
      
      fbtag(1) = "M2"
      fbfreq(1) = 0.000140518917083d0
      fbnfact(1) = 1d0
      fbeq(1) = 0d0            
      
      nfrbou = 0
      DO i = 1,nbou
        
        segtype = fbseg(2,i)
        nbnds = fbseg(1,i)          
          
        IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
          nfrbou = nfrbou + 1
          fbseg(1,nfrbou) = fbseg(1,i)
        ENDIF                    
      ENDDO      
      
      DO bou = 1,nfrbou
        nbnds = fbseg(1,bou)
        DO nd = 1,nbnds
          DO bfr = 1,nfbfr
            fbamp(nd,bou,bfr) = 1d0
            fbph(nd,bou,bfr) = 0d0
          ENDDO
        ENDDO        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      nfbsfr = 1
      CALL alloc_forcing_arrays(3)
      
      fbstag(1) = "SURGE"
      fbsbgn(1) = 1.5d0
      fbsend(1) = 3.0d0
      fbssig(1) = 0.18d0
      
      DO bou = 1,nfrbou
        nbnds = fbseg(1,bou)
        DO nd = 1,nbnds
          DO bfr = 1,nfbsfr
            fbsamp(nd,bou,bfr) = 16d0
          ENDDO
        ENDDO        
      ENDDO      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
 
      
      OPEN(UNIT=15, FILE=TRIM(ADJUSTL(bou_file)))
            
      WRITE(15,"(I7)") nobfr
      
      DO bfr = 1,nobfr
        WRITE(15,"(A10)") obtag(bfr)
        WRITE(15,"(3(E24.17,1x))") obfreq(bfr),obnfact(bfr),obeq(bfr)
      ENDDO
      
      DO bfr = 1,nobfr
        WRITE(15,"(A10)") obtag(bfr) 
        DO bou = 1,nope
          DO nd = 1,obseg(bou)
            WRITE(15,"(2(E24.17,1x))") obamp(nd,bou,bfr),obph(nd,bou,bfr)
          ENDDO
        ENDDO
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      WRITE(15,"(I7)") nfbfr
      
      DO bfr = 1,nfbfr
        WRITE(15,"(A10)") fbtag(bfr)
        WRITE(15,"(3(E24.17,1x))") fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
      ENDDO
      
      DO bfr = 1,nfbfr
        WRITE(15,"(A10)") fbtag(bfr) 
        DO bou = 1,nfrbou
          DO nd = 1,fbseg(1,bou)
            IF (mod(nd,10) == 0) THEN
              WRITE(15,"(1x,2(E24.17,1x))") fbamp(nd,bou,bfr),fbph(nd,bou,bfr)            
            ELSE
              WRITE(15,"(2(E24.17,1x))") fbamp(nd,bou,bfr),fbph(nd,bou,bfr)
            ENDIF
          ENDDO
        ENDDO
      ENDDO     
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      WRITE(15,"(I7)") nfbsfr
      
      DO bfr = 1,nfbsfr
        WRITE(15,"(A10)") fbstag(bfr)
        WRITE(15,"(3(E24.17,1x))") fbsbgn(bfr),fbsend(bfr),fbssig(bfr)
      ENDDO
      
      DO bfr = 1,nfbsfr
        WRITE(15,"(A10)") fbstag(bfr) 
        DO bou = 1,nfrbou
          DO nd = 1,fbseg(1,bou)
            WRITE(15,"(E24.17)") fbsamp(nd,bou,bfr)
          ENDDO
        ENDDO
      ENDDO          
      
      CLOSE(15)
      
      END PROGRAM write_bou_file