      SUBROUTINE modal2nodal()

      USE globals, ONLY: pres,nel_type,p,nnds,ndof,m2n,out_direc
      USE basis, ONLY: jacobi

      IMPLICIT NONE
      
      INTEGER :: i,j,m
      INTEGER :: et,pt
      INTEGER :: mndof,mnnds
      INTEGER :: nvert
      REAL(pres) :: r(4),s(4)
      REAL(pres) :: a(3),b(3)
      REAL(pres) :: Pi(4),Pj(4)      
      
      mndof = maxval(ndof)
      mnnds = maxval(nnds)
      
      ALLOCATE(m2n(mnnds,mndof,nel_type))      
      
      OPEN(unit=111,file=trim(out_direc) // "modal2nodal.d")
      
      DO et = 1,nel_type
        
        IF (mod(et,2) == 1) THEN
        
          nvert = 3
        
          r(1) = -1d0
          r(2) =  1d0
          r(3) = -1d0
          
          s(1) = -1d0
          s(2) = -1d0
          s(3) =  1d0
          
          DO pt = 1,nvert 
            IF (abs(s(pt) - 1d0) > 1d-14) THEN
              a(pt) = 2d0*(1d0+r(pt))/(1d0-s(pt))-1d0 
            ELSE 
              a(pt) = -1d0
            ENDIF
            b(pt) = s(pt)
          ENDDO       
          
          m = 0
          DO i = 0,p
            DO j = 0,p-i
            
              m = m + 1
              
              CALL jacobi(0    ,0,i,a,3,Pi)
              CALL jacobi(2*i+1,0,j,b,3,Pj)
              
              DO pt = 1,nvert
                m2n(pt,m,et) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
              ENDDO
              
            ENDDO
          ENDDO  
          

        ELSE IF (mod(et,2) == 0) THEN
        
          nvert = 4
        
          r(1) = -1d0
          r(2) =  1d0
          r(3) =  1d0
          r(4) = -1d0 
          
          s(1) = -1d0
          s(2) = -1d0
          s(3) =  1d0        
          s(4) =  1d0
          
          m = 0
          DO i = 0,p
            DO j = 0,p
            
              m = m + 1
              
              CALL jacobi(0,0,i,r,4,Pi)
              CALL jacobi(0,0,j,s,4,Pj)
              
              DO pt = 1,nvert 
                m2n(pt,m,et) = 2d0*Pi(pt)*Pj(pt)
              ENDDO
              
            ENDDO
          ENDDO
          
        ENDIF
        
        WRITE(111,"(I5,I5)") nvert, ndof(et)
        DO pt = 1,nvert
          WRITE(111,"(160(e24.17,1x))") (m2n(pt,m,et), m = 1,ndof(et))
        ENDDO        
          
      ENDDO
        
        
      

      RETURN
      END SUBROUTINE modal2nodal