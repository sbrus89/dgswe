      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "d008fbb5a186f02ee61f178479e2be916aa89957 +" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version