      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "eb9f136917f3dd44075f1caa2e4f1827a7412601 +" 
            
      PRINT*, "Version Information"
      PRINT*, "  Branch: ", gitBranch 
      PRINT*, "  SHA: ", gitSHA 
      PRINT*, " "
 
      END SUBROUTINE version