      MODULE transformation

      CONTAINS

      SUBROUTINE element_transformation(nnd,x,y,l,xpt,ypt,dldr,dlds,drdx,drdy,dsdx,dsdy,detJ)

      USE globals, ONLY: rp

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

      SUBROUTINE cpp_transformation(ypt,Sp)

      USE globals, ONLY: rp,r_earth
      USE read_dginp, ONLY: coord_sys,sphi0

      IMPLICIT NONE

      REAL(rp), INTENT(IN) :: ypt
      REAL(rp), INTENT(OUT) :: Sp

      IF (coord_sys == 1) THEN
        Sp = 1d0
      ELSE
        Sp = cos(sphi0)/cos(ypt/r_earth)        
      ENDIF      


      END SUBROUTINE cpp_transformation

      END MODULE transformation