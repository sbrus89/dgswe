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
      gitSHA = "9a71d760e96eb9f353fd80c4c6c5fdd59ca3ee9d +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-Iodir/ -O3 -g -traceback -C" 
      modified_files = "../src/plot_mod.F90 ../src/read_plot_input.F90" 
      compile_date = "Fri Mar  2 10:20:11 EST 2018" 
      host = "chl-tilos" 

      
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