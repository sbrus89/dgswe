      MODULE version
      
      IMPLICIT NONE
      
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitBranch = "master" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitSHA = "198b073bb862652bbd923489e7cc9bd9039731b1 +" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_version = "ifort version 14.0.0" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_flags = "-O3 -xHost -Iodir_dgprep/" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: modified_files = "../plot/plot_PE_grid.m ../plot/plot_grid.m ../plot/read_solution.m ../post/dgpost.F90 ../rimls/work/rimls.inp ../src/output.F90 ../src/read_input.F90 ../src/read_write_output.F90 ../src/system.F90 ../src/version.F90 ../stations/plot/plot_adcirc_stations.m ../stations/plot/plot_grid.m ../stations/plot/plot_stations.m ../stations/src/stations.f90 ../util/Makefile ../util/adjust Makefile dgswe.inp" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compile_date = "Thu Jun 16 16:40:18 EDT 2016" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: host = "chl-tilos" 
            
      CONTAINS      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE version_information(unit)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: unit
      INTEGER :: length,ind,sind,eind            

      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
            
      CALL directed_output("Version Information",unit)
      CALL directed_output("  Branch: "//TRIM(ADJUSTL(gitBranch)),unit)
      CALL directed_output("  SHA: "//TRIM(ADJUSTL(gitSHA)),unit)
      CALL directed_output("  Modified files: ",unit)

        
      sind = 1        
      ind = INDEX(modified_files," ")        
      eind = ind
      DO WHILE (ind > 0)
        CALL directed_output("    " //TRIM(ADJUSTL(modified_files(sind:eind))),unit)
        sind = eind+1
        ind = INDEX(modified_files(sind:length)," ")
        eind = eind + ind
      ENDDO
        
      CALL directed_output(" ",unit)
      CALL directed_output("Compiler Information",unit)
      CALL directed_output("  Version: "//TRIM(ADJUSTL(compiler_version)),unit)
      CALL directed_output("  Flags: "//TRIM(ADJUSTL(compiler_flags)),unit)
      CALL directed_output("  Date compiled: "//TRIM(ADJUSTL(compile_date)),unit)
      CALL directed_output("  Host: "//TRIM(ADJUSTL(host)),unit)
      CALL directed_output(" ",unit)
     
 
      END SUBROUTINE version_information
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     

      
      END MODULE version