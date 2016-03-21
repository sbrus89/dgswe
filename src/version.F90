      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      
      gitBranch = "master" 
      gitSHA = "9c70d5f9ac2be2f926eb5632c28a7eb503ba5c08 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-132 -heap-arrays -Iodir/ -traceback -g -C" 
      modified_files = "../src/evaluate.f90 ../src/find_element.f90 ../src/read_grid.f90 ../src/rimls.f90 Makefile" 
      compile_date = "Mon Mar 21 15:49:35 EDT 2016" 
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