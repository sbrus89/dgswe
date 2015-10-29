      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "104e1b957bd07605ab861eff5a83b334a099ef07 +" 
      compiler_version = "mpifort for MPICH version 3.1.4 ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -DCMPI -Iodir_dgswe_mpi/" 
      modified_files = "../grids/dble_prec_grid.F90 ../plot/plot_sol_fill.m ../src/dgswe.F90 ../src/edge_integration.F90 ../src/globals.F90 ../src/interp_forcing.F90 ../src/messenger2.F90 ../src/read_grid.F90 ../src/version.F90 Makefile dgswe.inp" 
      compile_date = "Thu Oct 29 10:50:30 EDT 2015" 
      host = "chl-tilos" 
      
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