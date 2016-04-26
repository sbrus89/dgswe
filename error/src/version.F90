      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "master" 
      gitSHA = "5e75f1a147deeb4ad4367f70ace31d06e25a5ea9 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -132 -Iodir/ " 
      modified_files = "error.inp" 
      compile_date = "Fri Nov  6 14:40:40 EST 2015" 
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