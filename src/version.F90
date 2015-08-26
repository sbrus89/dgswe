      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "8429a7d25b85b789be4b313d8e2599e8b9595a44 +" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version