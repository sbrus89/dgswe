      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "master" 
      gitSHA = "b65c83c5ab2f0ad91855295bd2392a3dcd3d5639 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-132 -C -g -traceback -Iodir/ " 
      modified_files = "../src/allocation.f90 ../src/bathy_interp.F90 ../src/evaluate.F90 ../src/find_element.f90 ../src/globals.f90 ../src/read_grid.f90 ../src/version.F90 Makefile bathy.inp" 
      compile_date = "Mon Dec  7 15:39:57 EST 2015" 
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