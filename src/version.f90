      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      
      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "master" 
      gitSHA = "d22a699afc37a05d152809d5ccc32a561a82b648 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-132 -heap-arrays -Iodir/ " 
      modified_files = "../src/evaluate.f90 ../src/version.f90" 
      compile_date = "Sun Nov  8 15:54:10 EST 2015" 
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