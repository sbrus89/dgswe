      SUBROUTINE linear_solve(et,sel,eel,ndof)

      USE globals, ONLY: rhsZ,rhsQx,rhsQy, &
                         MirhsZ,MirhsQx,MirhsQy, &
                         mmi

      IMPLICIT NONE

      INTEGER :: et,sel,eel,ndof
      INTEGER :: i,j,l,el,m

            SELECT CASE(et)
               
              CASE(1)
                DO l = 1,ndof
                  DO el = sel,eel
                    MirhsZ(el,l)  = MirhsZ(el,l)  + mmi(el,1)*rhsZ(el,l) 
                    MirhsQx(el,l) = MirhsQx(el,l) + mmi(el,1)*rhsQx(el,l) 
                    MirhsQy(el,l) = MirhsQy(el,l) + mmi(el,1)*rhsQy(el,l) 
                    
                    rhsZ(el,l) = 0d0
                    rhsQx(el,l) = 0d0
                    rhsQy(el,l) = 0d0
                  ENDDO
                ENDDO
 
              CASE DEFAULT
              
                m = 1
                DO i = 1,ndof
                  DO j = 1,ndof
                    DO el = sel,eel
                      MirhsZ(el,i)  = MirhsZ(el,i)  + mmi(el,m)*rhsZ(el,j) 
                      MirhsQx(el,i) = MirhsQx(el,i) + mmi(el,m)*rhsQx(el,j) 
                      MirhsQy(el,i) = MirhsQy(el,i) + mmi(el,m)*rhsQy(el,j)
                      
                      rhsZ(el,j) = 0d0
                      rhsQx(el,j) = 0d0
                      rhsQy(el,j) = 0d0                      
                    ENDDO
                    m = m + 1
                  ENDDO
                ENDDO            
            END SELECT      

      END SUBROUTINE linear_solve