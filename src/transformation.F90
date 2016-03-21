      MODULE transformation
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      

      CONTAINS
      
      SUBROUTINE init_element_coordinates(ne,mnnds,el_type,nverts,xy,ect,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: mnnds
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: elxy
            
      INTEGER :: el,nd
      INTEGER :: et,nv
      INTEGER :: alloc_status
      
      ALLOCATE(elxy(mnnds,ne,2), STAT=alloc_status)
      
      
      DO el = 1,ne        
        et = el_type(el)
        nv = nverts(et)
        DO nd = 1,nv
          elxy(nd,el,1) = xy(1,ect(nd,el))
          elxy(nd,el,2) = xy(2,ect(nd,el))
        ENDDO              
      ENDDO          
      
      RETURN
      END SUBROUTINE init_element_coordinates
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE element_transformation(nnd,x,y,l,xpt,ypt,dldr,dlds,drdx,drdy,dsdx,dsdy,detJ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nnd

      REAL(rp), DIMENSION(:), INTENT(IN) :: x,y
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: xpt,ypt
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(OUT), OPTIONAL :: drdx,drdy,dsdx,dsdy,detJ

      INTEGER :: nd
      REAL(rp) :: dxdr,dxds,dydr,dyds
      INTEGER :: calc_inv


      IF ( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
           PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
           PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. &
           PRESENT(detJ) ) THEN
        calc_inv = 1
      ELSE
        calc_inv = 0
      ENDIF
     

      xpt = 0d0
      ypt = 0d0
          
      DO nd = 1,nnd        
              
        xpt = xpt + l(nd)*x(nd)
        ypt = ypt + l(nd)*y(nd)
                      
      ENDDO

         
      IF (calc_inv == 1) THEN

        dxdr = 0d0
        dxds = 0d0
        dydr = 0d0
        dyds = 0d0
            
          
        DO nd = 1,nnd
        
          dxdr = dxdr + dldr(nd)*x(nd)
          dxds = dxds + dlds(nd)*x(nd)
          dydr = dydr + dldr(nd)*y(nd)
          dyds = dyds + dlds(nd)*y(nd)              
                      
        ENDDO
            
        detJ = dxdr*dyds - dxds*dydr
            
        drdx =  dyds/detJ
        drdy = -dxds/detJ
        dsdx = -dydr/detJ
        dsdy =  dxdr/detJ

      ENDIF

      END SUBROUTINE element_transformation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     


      END MODULE transformation