      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "e8aa51ab1266615ae0bb08d35bc6adbaaa513255 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -Iodir_dgpost/" 
      modified_files = "../grids/converge.grd ../grids/dble_prec_grid.F90 ../grids/refine_grid.m ../plot/plot_sol_fill.m ../post/dgpost.F90 ../src/version.F90 Makefile" 
      compile_date = "Thu Nov 19 13:41:44 EST 2015" 
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