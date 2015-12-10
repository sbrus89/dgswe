      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "master" 
      gitSHA = "80b84f32c43c755f2fd91910ca151b8f06f07bef +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-132 -C -g -traceback -Iodir/ " 
      modified_files = "../src/read_input.f90 ../src/version.F90 ../src/write_results.f90 bathy.inp" 
      compile_date = "Thu Dec 10 17:54:38 EST 2015" 
      host = "chl-tilos" 
      
      length = LEN(TRIM(ADJUSTL(modified_files)))  
            
        PRINT*, ""    
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