      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "022fffcb5a6cbee435671f56582a9a10a81adbc0 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -Iodir_dgswe/" 
      modified_files = "../src/globals.F90 ../src/version.F90 Makefile" 
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", TRIM(ADJUSTL(gitBranch))
        PRINT*, "  SHA: ", TRIM(ADJUSTL(gitSHA)) 
        PRINT*, "  Modified files: ", TRIM(ADJUSTL(modified_files))
        PRINT*, " "
        PRINT*, "Compiler Information"
        PRINT*, "  Version: ", TRIM(ADJUSTL(compiler_version))
        PRINT*, "  Flags: ", TRIM(ADJUSTL(compiler_flags))
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version