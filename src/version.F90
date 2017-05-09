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
      gitSHA = "3e4fc47e7c95982d51e7002611a11f7165bcca6e +" 
      compiler_version = "ifort version 14.0.0" 
      compiler_flags = "-O3 -xHost -Iodir_dgswe/" 
      modified_files = "" 
      compile_date = "Tue May  9 15:07:59 EDT 2017" 
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