      PROGRAM dgpost
      
      USE globalize, ONLY: globalize_solution

      IMPLICIT NONE
      
      INTEGER :: i,pe,el,dof
      INTEGER :: npe,tne,mne,ne,ndof,lines
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER, PARAMETER :: pres = kind(1d0)
      INTEGER :: lname
      CHARACTER(6) :: dirname
      CHARACTER(100) :: file_name
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: H_global, Qx_global, Qy_global
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: hb_local, hb_global
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: t
      
      CALL version()

      OPEN(UNIT=19,FILE='PE0000/fort.80')      
      READ(19,*) npe
      CLOSE(19)

      PRINT*, "Globalizing solutions"
      DO pe = 1,npe
      
        dirname = "PE0000"
        lname = 6
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1
        
        OPEN(UNIT=19,FILE=dirname(1:lname)//'/'//'fort.80',POSITION='rewind')
        
        READ(19,*)
        READ(19,*) tne
        READ(19,*) mne
        READ(19,*) ne
        READ(19,*) ndof
        READ(19,*) lines
        
        IF(.not. ALLOCATED(lel2gel)) THEN
          ALLOCATE(lel2gel(mne))
        ENDIF
        
        DO el = 1,ne
          READ(19,*) i, lel2gel(el)
        ENDDO
        
        CLOSE(19)   
        
        
        
        IF(.not. ALLOCATED(hb_local)) THEN
          ALLOCATE(hb_local(mne,ndof))
        ENDIF        
        
        IF(.not. ALLOCATED(H_global))THEN
          ALLOCATE(t(lines+2))
          ALLOCATE(H_global(tne,ndof,lines+1),Qx_global(tne,ndof,lines+1),Qy_global(tne,ndof,lines+1),hb_global(tne,ndof))
        ENDIF  
        
        file_name = dirname(1:lname)//'/'//'solution_H.d'
        CALL globalize_solution(file_name,ne,mne,ndof,lel2gel,t,H_global)
        
        file_name = dirname(1:lname)//'/'//'solution_Qx.d'
        CALL globalize_solution(file_name,ne,mne,ndof,lel2gel,t,Qx_global)   
        
        file_name = dirname(1:lname)//'/'//'solution_Qy.d'
        CALL globalize_solution(file_name,ne,mne,ndof,lel2gel,t,Qy_global)        



        OPEN(UNIT=65,FILE=dirname(1:lname)//'/'//'hb_modal.d')
        DO dof = 1,4
          READ(65,6364) (hb_local(el,dof), el = 1,ne)
        ENDDO
          
        DO dof = 1,4
          DO el = 1,ne
            hb_global(lel2gel(el),dof) = hb_local(el,dof)
          ENDDO
        ENDDO
          
        CLOSE(642)        

        
      ENDDO
      PRINT*, "  done"      
      
      
      
      PRINT*, "Writing solutions"
      OPEN(UNIT=63,FILE='solution_H.d')
      WRITE(63,*) 'H solution'
      DO tstep = 1,lines+1
        WRITE(63,*) t(tstep)
        DO dof = 1,ndof
          WRITE(63,6364) (H_global(el,dof,tstep), el=1,tne)
        ENDDO
      ENDDO
      CLOSE(63)
      
      OPEN(UNIT=641,FILE='solution_Qx.d')
      WRITE(641,*) 'Qx solution'
      DO tstep = 1,lines+1
        WRITE(641,*) t(tstep)
        DO dof = 1,ndof
          WRITE(641,6364) (Qx_global(el,dof,tstep), el=1,tne)
        ENDDO
      ENDDO
      CLOSE(641)
      
      OPEN(UNIT=642,FILE='solution_Qy.d')
      WRITE(642,*) 'Qy solution'
      DO tstep = 1,lines+1
        WRITE(642,*) t(tstep)
        DO dof = 1,ndof
          WRITE(642,6364) (Qy_global(el,dof,tstep), el=1,tne)
        ENDDO
      ENDDO
      CLOSE(642)
      
      OPEN(UNIT=65,FILE='hb_modal.d')
      DO dof = 1,ndof
        WRITE(65,6364) (hb_global(el,dof), el=1,tne)
      ENDDO
      CLOSE(65)
      
      PRINT*, "  done"
      
      CALL system('cp PE0000/modal2nodal.d .')
      CALL system('cp PE0000/projection.d .')
      
 6364  FORMAT(160000(e24.17,1x))
      END PROGRAM dgpost
