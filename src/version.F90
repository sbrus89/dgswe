      SUBROUTINE version()
      
      USE globals, ONLY: gitSHA,gitBranch,compiler_version,compiler_flags,modified_files,compile_date,host

      IMPLICIT NONE
      
      INTEGER :: length,ind,sind,eind
      
      
      gitBranch = "master" 
      gitSHA = "87e1b75aecedf62fec85c296af422732b6dd4c38 +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -Iodir_dgprep/" 
      modified_files = "../bathy_interp/work/bathy.inp ../error/work/error.inp ../grids/convert14.py ../grids/convert2dm.py ../grids/dble_prec_grid.F90 ../grids/refine_grid.m ../plot/drawNGonMesh4.m ../plot/plot_PE_grid.m ../plot/plot_grid.m ../plot/plot_sol_fill.m ../rimls/plot/drawNGonMesh4.m ../rimls/plot/plot_nodes.m ../rimls/src/write_results.f90 ../rimls/work/rimls.inp ../spline/plot/plot_spline_constant.m ../spline/src/calc_spline.f90 ../spline/src/spline.f90 ../spline/work/Makefile ../spline/work/spline.inp ../src/bathymetry_interp_mod.F90 ../src/grid_file_mod.F90 ../src/version.F90 ../stations/src/stations.f90 ../stations/work/Makefile dgswe.inp" 
      compile_date = "Wed May  4 02:51:15 EDT 2016" 
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