      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "10d6ff9cf6a2569f92ee73a9384393d292f46d38 +" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version