      MODULE version
      
      IMPLICIT NONE
      
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitBranch = "master" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitSHA = "b08fa604cd455d8dd1e05739072791b1f8c07833 +" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_version = "ifort version 14.0.0" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_flags = "-O3 -xHost -Iodir_dgprep/" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: modified_files = "../plot/plot_PE_grid.m ../plot/plot_grid.m ../rimls/work/rimls.inp ../stations/plot/plot_adcirc_stations.m ../stations/plot/plot_grid.m ../stations/plot/plot_stations.m ../stations/src/stations.f90" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compile_date = "Mon Jun 13 14:23:46 EDT 2016" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: host = "chl-tilos" 
            
      CONTAINS      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE version_information(unit)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: unit
      INTEGER :: length,ind,sind,eind            

      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
            
      CALL output("Version Information",unit)
      CALL output("  Branch: "//TRIM(ADJUSTL(gitBranch)),unit)
      CALL output("  SHA: "//TRIM(ADJUSTL(gitSHA)),unit)
      CALL output("  Modified files: ",unit)

        
      sind = 1        
      ind = INDEX(modified_files," ")        
      eind = ind
      DO WHILE (ind > 0)
        CALL output("    " //TRIM(ADJUSTL(modified_files(sind:eind))),unit)
        sind = eind+1
        ind = INDEX(modified_files(sind:length)," ")
        eind = eind + ind
      ENDDO
        
      CALL output(" ",unit)
      CALL output("Compiler Information",unit)
      CALL output("  Version: "//TRIM(ADJUSTL(compiler_version)),unit)
      CALL output("  Flags: "//TRIM(ADJUSTL(compiler_flags)),unit)
      CALL output("  Date compiled: "//TRIM(ADJUSTL(compile_date)),unit)
      CALL output("  Host: "//TRIM(ADJUSTL(host)),unit)
      CALL output(" ",unit)
     
 
      END SUBROUTINE version_information
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
      SUBROUTINE output(string,unit)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: unit
      
      IF (unit == 6) THEN
        PRINT*, string
      ELSE
        WRITE(unit,*) string
      ENDIF
      
      
      
      RETURN
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      END MODULE version