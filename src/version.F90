      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      
      gitBranch = "master" 
      gitSHA = "3ac8d418c4e072e6befcb797a9f4a1a2e6e41f45 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-C -g -traceback -Iodir_dgswe/" 
      modified_files = "../src/curvilinear.F90 ../src/element_data.F90 ../src/grid_file_mod.F90 ../src/version.F90 Makefile" 
      compile_date = "Sun Mar 20 20:04:19 EDT 2016" 
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