      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host
      USE messenger2, ONLY: myrank

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      gitBranch = "quadtri_mixed_mpi" 
      gitSHA = "b9148cceb570a4ef5327891f54f994508b6ba2f9 +" 
      compiler_version = "mpifort for MPICH version 3.1.4 ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -DCMPI -Iodir_dgswe_mpi/" 
      modified_files = "../grids/dble_prec_grid.F90 ../plot/plot_sol_fill.m ../post/dgpost.F90 ../prep/src/decomp2.F90 ../prep/src/write_files.F90 ../src/curvilinear.F90 ../src/version.F90 Makefile dgswe.inp" 
      compile_date = "Mon Nov  2 17:39:32 EST 2015" 
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