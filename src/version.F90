      MODULE version
      
      IMPLICIT NONE
      
      CHARACTER(:), ALLOCATABLE :: gitBranch
      CHARACTER(:), ALLOCATABLE :: gitSHA 
      CHARACTER(:), ALLOCATABLE :: compiler_version 
      CHARACTER(:), ALLOCATABLE :: compiler_flags
      CHARACTER(:), ALLOCATABLE :: modified_files
      CHARACTER(:), ALLOCATABLE :: compile_date
      CHARACTER(:), ALLOCATABLE :: host
            
      CONTAINS      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE version_information(unit)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: unit
      INTEGER :: length,ind,sind,eind        

      gitBranch = "master" 
      gitSHA = "b6b7d91a77e0e05c32864f555288a851847fe10d +" 
      compiler_version = "ifort version 12.1.3" 
      compiler_flags = "-Iodir/ -g -traceback -C" 
      modified_files = "../src/evaluate_mod.F90 ../src/initialize.F90 ../src/plot_dg.F90 ../src/plot_globals.F90 ../src/plot_mod.F90 Makefile" 
      compile_date = "Thu Dec  7 16:10:13 EST 2017" 
      host = "sbrus-laptop" 

      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
            
      CALL directed_output("Version Information",unit)
      CALL directed_output("  Branch: "//TRIM(ADJUSTL(gitBranch)),unit)
      CALL directed_output("  SHA: "//TRIM(ADJUSTL(gitSHA)),unit)
!       CALL directed_output("  Modified files: ",unit)
! 
!         
!       sind = 1        
!       ind = INDEX(modified_files," ")        
!       eind = ind
!       DO WHILE (ind > 0)
!         CALL directed_output("    " //TRIM(ADJUSTL(modified_files(sind:eind))),unit)
!         sind = eind+1
!         ind = INDEX(modified_files(sind:length)," ")
!         eind = eind + ind
!       ENDDO
        
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