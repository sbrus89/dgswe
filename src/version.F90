      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "ef2a23106e648a13c9d864c8c68c9f7fc713d9a4 +" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version