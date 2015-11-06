      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "e718172ce2f3006c1d8c4a4db48116ae2a023a78 +" 
      compiler_version = "ifort version 15.0.4" 
      compiler_flags = "-O3 -xHost -Iodir_dgprep/" 
      modified_files = "" 
      compile_date = "Tue Nov  3 14:53:21 EST 2015" 
      host = "aegaeon crc nd edu" 
      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
      IF (myrank == 0) THEN     
            
        PRINT*, "Version Information"
        PRINT*, "  Branch: ", TRIM(ADJUSTL(gitBranch))
        PRINT*, "  SHA: ", TRIM(ADJUSTL(gitSHA)) 
        PRINT*, "  Modified files: "

        
        sind = 1        
        ind = INDEX(modified_files," ")        
        eind = ind
        DO WHILE (ind > 0)
          PRINT*, "    " ,TRIM(ADJUSTL(modified_files(sind:eind)))
          sind = eind+1
          ind = INDEX(modified_files(sind:length)," ")
          eind = eind + ind
        ENDDO
        
        PRINT*, " "
        PRINT*, "Compiler Information"
        PRINT*, "  Version: ", TRIM(ADJUSTL(compiler_version))
        PRINT*, "  Flags: ", TRIM(ADJUSTL(compiler_flags))
        PRINT*, "  Date compiled: ", TRIM(ADJUSTL(compile_date))
        PRINT*, "  Host: ", TRIM(ADJUSTL(host))
        PRINT*, " "
      
      ENDIF
 
      END SUBROUTINE version