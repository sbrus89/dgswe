      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "c1f188714ef321d0f44a1eaf980b6da0f271e05a +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O2 -align array32byte -align rec32byte -Iodir_dgswe/" 
      
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", gitBranch 
        PRINT*, "  SHA: ", gitSHA 
        PRINT*, " "
        PRINT*, "Compiler Information"
        PRINT*, "  Version: ", compiler_version
        PRINT*, "  Flags: ", compiler_flags
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version