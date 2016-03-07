      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      
      gitBranch = "master" 
      gitSHA = "7a581de70b24588143cc802e618021312cf97e2a +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-C -g -traceback -Iodir_dgswe/" 
      modified_files = "../rimls/src/version.f90 ../src/basis.F90 ../src/modal2nodal.F90 ../src/version.F90 Makefile" 
      compile_date = "Mon Mar  7 16:47:36 EST 2016" 
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