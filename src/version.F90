      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "a90514ca26982fbef203aab77c5aa0913acaa9c3 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -Iodir_dgswe/" 
      
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", TRIM(ADJUSTL(gitBranch))
        PRINT*, "  SHA: ", TRIM(ADJUSTL(gitSHA)) 
        PRINT*, " "
        PRINT*, "Compiler Information"
        PRINT*, "  Version: ", TRIM(ADJUSTL(compiler_version))
        PRINT*, "  Flags: ", TRIM(ADJUSTL(compiler_flags))
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version