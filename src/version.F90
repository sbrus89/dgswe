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
      gitSHA = "8b555ec4b6c5c54c6b8b3c4802a7908430a17bbc +"
      compiler_version = "Using built-in specs. COLLECT_GCC=gfortran"
      compiler_flags =  "-O3  -Iodir_dgswe/ -DMAC"
      modified_files = "../src/plot_mod.F90 ../src/read_plot_input.F90" 
      compile_date =  "Wed Jun  6 15:32:23 MDT 2018"
      host =  "pn1808279 lanl gov"

      
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
