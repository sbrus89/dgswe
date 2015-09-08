      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "5ec7f543a31f042f556688d5a19e404523f28bea +" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version