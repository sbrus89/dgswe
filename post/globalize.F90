      MODULE globalize
      
      USE globals, ONLY:rp
      
      CONTAINS
      
      SUBROUTINE globalize_solution(file_name,ne,mne,ndof,lel2gel,t,sol_global)

      IMPLICIT NONE
      
      CHARACTER(100), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: mne
      INTEGER, INTENT(IN) :: ndof
      INTEGER, DIMENSION(:) :: lel2gel
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: sol_global
      
      INTEGER :: el,dof
      INTEGER :: tstep
      INTEGER :: read_stat
      REAL(rp) :: sol_local(mne,ndof)

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

      END MODULE globalize