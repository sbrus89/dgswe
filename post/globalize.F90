      MODULE globalize
      
      USE globals, ONLY:rp
      
      CONTAINS
      
      SUBROUTINE globalize_solution(file_name,ne,ndof,lel2gel,t,sol_global)

      IMPLICIT NONE
      
      CHARACTER(100), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: ndof
      INTEGER, DIMENSION(:) :: lel2gel
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: sol_global
      
      INTEGER :: el,dof
      INTEGER :: tstep
      INTEGER :: read_stat
      REAL(rp) :: sol_local(ne,ndof)

      OPEN(UNIT=63,FILE=TRIM(file_name))
      READ(63,*)
      tstep = 0
      DO         
        READ(63,*,IOSTAT=read_stat) t(tstep + 1)
        IF(read_stat < 0) THEN
          EXIT 
        ENDIF          
  
        DO dof = 1,ndof
          READ(63,6364) (sol_local(el,dof), el = 1,ne)   
        ENDDO
          
        tstep = tstep + 1
        
        DO dof = 1,ndof
          DO el = 1,ne
            sol_global(lel2gel(el),dof,tstep) = sol_local(el,dof)
          ENDDO
        ENDDO
                
      ENDDO 
      CLOSE(63)

 6364  FORMAT(160000(e24.17,1x))

      RETURN
      END SUBROUTINE globalize_solution
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      SUBROUTINE globalize_stations(file_name,nsta,sta_l2g,t,sta_global)

      IMPLICIT NONE
      
      CHARACTER(100), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: nsta
      INTEGER, DIMENSION(:) :: sta_l2g
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: sta_global
      
      INTEGER :: sta
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: nsta2
      REAL(rp) :: sta_local(nsta)

      OPEN(UNIT=61,FILE=TRIM(file_name))
      READ(61,*)
      READ(61,*) nsta2
      
      IF (nsta /= nsta2) THEN
        PRINT*, "Number of local stations do not match"
        STOP
      ENDIF
      
      tstep = 0
      DO         
        READ(61,*,IOSTAT=read_stat) t(tstep + 1)
        IF(read_stat < 0) THEN
          EXIT 
        ENDIF          
  
        DO sta = 1,nsta
          READ(61,*) sta_local(sta)   
        ENDDO
          
        tstep = tstep + 1
        
          DO sta = 1,nsta
            sta_global(sta_l2g(sta),tstep) = sta_local(sta)
          ENDDO
              
      ENDDO 
      CLOSE(61)


      RETURN
      END SUBROUTINE globalize_stations
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      SUBROUTINE write_global_sol(file_name,nout_sol,ne,ndof,t_sol,sol_global)
      
      IMPLICIT NONE
      
      CHARACTER(100), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: nout_sol
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: ndof
      REAL(rp), DIMENSION(:), INTENT(IN) :: t_sol
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: sol_global
      
      INTEGER :: tstep
      INTEGER :: dof
      INTEGER :: el
      
      
      OPEN(UNIT=63,FILE=TRIM(file_name))
      WRITE(63,*) 'Globalized Solution'
      DO tstep = 1,nout_sol+1
        WRITE(63,*) t_sol(tstep)
        DO dof = 1,ndof
          WRITE(63,6364) (sol_global(el,dof,tstep), el=1,ne)
        ENDDO
      ENDDO
      CLOSE(63)      
      
 6364  FORMAT(160000(e24.17,1x))      
      
      RETURN
      END SUBROUTINE write_global_sol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_global_sta(file_name,nout_sta,nsta,t_sta,sta_global)
      
      IMPLICIT NONE
      
      CHARACTER(100), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: nout_sta
      INTEGER, INTENT(IN) :: nsta
      REAL(rp), DIMENSION(:), INTENT(IN) :: t_sta
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sta_global
      
      INTEGER :: tstep
      INTEGER :: sta
      
      OPEN(UNIT=61,FILE=TRIM(file_name))
      WRITE(61,*) 'Globalized Stations'
      WRITE(61,*) nsta
      DO tstep = 1,nout_sta+1
        WRITE(61,*) t_sta(tstep)
        DO sta = 1,nsta
          WRITE(61,"(e24.17)") sta_global(sta,tstep)
        ENDDO
      ENDDO
      CLOSE(61)
          
      
      RETURN
      END SUBROUTINE write_global_sta
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE globalize