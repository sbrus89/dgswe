      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      
      gitBranch = "master" 
      gitSHA = "7a76f32a1143cc07247e1bf34b5b0e2273ce7a58 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-132 -heap-arrays -Iodir/ " 
      modified_files = "../src/connect.f90 ../src/evaluate.f90 ../src/globals.f90 ../src/version.f90 ../src/write_results.f90 Makefile rimls.inp" 
      compile_date = "Sun Mar 20 19:27:29 EDT 2016" 
      host = "chl-tilos" 
      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
            
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
     
 
      END SUBROUTINE version