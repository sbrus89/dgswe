      MODULE bathymetry_interp_mod

      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE bathymetry_interp_eval(n,hb,l,h,dldr,dlds,drdx,drdy,dsdx,dsdy,Sp,dhdx,dhdy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(IN) :: hb
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: h      
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(IN), OPTIONAL :: drdx,drdy,dsdx,dsdy
      REAL(rp), INTENT(IN), OPTIONAL :: Sp
      REAL(rp), INTENT(OUT), OPTIONAL :: dhdx,dhdy
      
      INTEGER :: nd
      INTEGER :: calc_deriv
      
      IF( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
          PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
          PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. & 
          PRESENT(Sp)   .AND.                     &
          PRESENT(dhdx) .AND. PRESENT(dhdy) ) THEN
          
        calc_deriv = 1
          
      ELSE
      
        calc_deriv = 0
          
      ENDIF

      h = 0d0      
      DO nd = 1,n                
        h =  h + l(nd)*hb(nd)                           
      ENDDO         
      
      IF (calc_deriv == 1) THEN
      
        dhdx = 0d0
        dhdy = 0d0
        DO nd = 1,n  
          dhdx = dhdx + (dldr(nd)*drdx + dlds(nd)*dsdx)*hb(nd)*Sp
          dhdy = dhdy + (dldr(nd)*drdy + dlds(nd)*dsdy)*hb(nd)              
        ENDDO      
      
      ENDIF 
      
      
      RETURN
      END SUBROUTINE bathymetry_interp_eval
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      END MODULE bathymetry_interp_mod