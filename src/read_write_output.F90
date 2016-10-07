      MODULE read_write_output
      
      USE globals, ONLY:rp      

      INTEGER :: unit_counter = 10
      
      CONTAINS     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE time_snaps(opt,snap,tf,dt,tskp,nout)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: opt
      REAL(rp), INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: tf
      REAL(rp), INTENT(IN) :: dt
      INTEGER, INTENT(OUT) :: tskp
      INTEGER, INTENT(OUT) :: nout
      
      
        ! Set number of timesteps between output
        SELECT CASE (opt)
          CASE (0)                    ! output off
            tskp = int(tf/dt) + 1     ! more than the number of total iterations
            nout = 0
          CASE (1)                    ! snap = number of total snaps
            tskp = int(tf/(snap*dt))  
            nout = int(snap)
          CASE (2)                    ! snap = number of timesteps between snaps
            tskp = int(snap)          
            nout = int(tf/(snap*dt))
          CASE (3)                    ! snap = number of seconds between snaps
            tskp = int(snap/dt)       
            nout = int(tf/snap)
          CASE DEFAULT
            PRINT*, "output option not supported"
        END SELECT          
      
      RETURN
      END SUBROUTINE time_snaps      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      SUBROUTINE file_init(out_direc,filename,nrow,ncol,nsnap,file_unit,post)
      
      
      USE version, ONLY: version_information
      USE read_dginp, ONLY: write_input,write_file_SHAs      
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nrow
      INTEGER, INTENT(IN) :: ncol
      INTEGER, INTENT(IN) :: nsnap    
      INTEGER, INTENT(OUT) :: file_unit  
      LOGICAL, INTENT(IN), OPTIONAL :: post
      
      INTEGER :: copy
      
      copy = 0
      IF (PRESENT(post)) THEN ! if post, copy version information from PE0000
        copy = 1
      ENDIF
               
      file_unit = unit_counter               
      OPEN(UNIT=file_unit,file=TRIM(ADJUSTL(out_direc)) // TRIM(ADJUSTL(filename)))
      unit_counter = unit_counter + 1
      
      IF (copy == 0) THEN
        CALL version_information(file_unit)
      ELSE
        CALL copy_version_header(file_unit,filename) 
      ENDIF
      WRITE(file_unit,"(A)") "-----------------------------------------------------------------------"      
      CALL write_input(file_unit)
      CALL write_file_SHAs(file_unit,out_direc)
      WRITE(file_unit,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      WRITE(file_unit,*) nrow,ncol,nsnap      

      RETURN
      END SUBROUTINE file_init   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      SUBROUTINE write_solution_snap(file_unit,nrow,ncol,trnspse,t_sol,solution)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: nrow      
      INTEGER, INTENT(IN) :: ncol
      CHARACTER(1), INTENT(IN) :: trnspse      
      REAL(rp), INTENT(IN) :: t_sol
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: solution
      
      INTEGER :: row
      INTEGER :: col
      INTEGER :: flip

      
      WRITE(file_unit,"(e24.17)") t_sol
      IF (trnspse == "T") THEN
        DO row = 1,nrow
          WRITE(file_unit,"(2x,*(e24.17,1x))") (solution(col,row), col=1,ncol)
        ENDDO
      ELSE
        DO row = 1,nrow
          WRITE(file_unit,"(2x,*(e24.17,1x))") (solution(row,col), col=1,ncol)
        ENDDO      
      ENDIF
      
      
      RETURN
      END SUBROUTINE write_solution_snap
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE write_solution_full(out_direc,file_name,ndof,ne,nsnap,trnspse,t_sol,solution,post)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: ndof
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nsnap
      CHARACTER(1), INTENT(IN) :: trnspse      
      REAL(rp), DIMENSION(:), INTENT(IN) :: t_sol
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: solution
      LOGICAL, INTENT(IN), OPTIONAL :: post
      
      INTEGER :: file_unit
      INTEGER :: snap
      
      
      IF (PRESENT(post)) THEN ! if post, copy version information from PE0000
        CALL file_init(out_direc,file_name,ndof,ne,nsnap,file_unit,post) 
      ELSE 
        CALL file_init(out_direc,file_name,ndof,ne,nsnap,file_unit) 
      ENDIF
            
     
      
      DO snap = 1,nsnap
        CALL write_solution_snap(file_unit,ndof,ne,trnspse,t_sol(snap),solution(:,:,snap))
      ENDDO
      
      CLOSE(file_unit)
      
      RETURN
      END SUBROUTINE write_solution_full      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE read_solution_header(out_direc,filename,nrow,ncol,nsnap,file_unit)
      
      IMPLICIT NONE
            
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(OUT) :: nrow
      INTEGER, INTENT(OUT) :: ncol
      INTEGER, INTENT(OUT) :: nsnap    
      INTEGER, INTENT(OUT) :: file_unit  
      
      CHARACTER(10) :: flag
      LOGICAL :: file_exists       
      
      INQUIRE(FILE=TRIM(ADJUSTL(out_direc)) // TRIM(ADJUSTL(filename)), EXIST=file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "file does not exist"
        STOP        
      ENDIF      
                           
      file_unit = unit_counter                               
      OPEN(UNIT=file_unit,file=TRIM(ADJUSTL(out_direc)) // TRIM(ADJUSTL(filename)))
      unit_counter = unit_counter + 1
      
  pre:DO 
      
        READ(file_unit,*) flag
        
        IF (flag == "!!!!!!!!!!") THEN
          EXIT pre
        ENDIF
      
      ENDDO pre
      
      READ(file_unit,*) nrow,ncol,nsnap        
      
      RETURN
      END SUBROUTINE read_solution_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE read_solution_snap(file_unit,nrow,ncol,trnspse,map,read_stat,t,solution)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit      
      INTEGER, INTENT(OUT) :: read_stat        
      INTEGER, INTENT(IN) :: nrow
      INTEGER, INTENT(IN) :: ncol      
      CHARACTER(1), INTENT(IN) :: trnspse      
      INTEGER, DIMENSION(:), INTENT(IN) :: map      
      REAL(rp), INTENT(OUT) :: t
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: solution   
      
      INTEGER :: col
      INTEGER :: row

      READ(file_unit,*,IOSTAT=read_stat) t
      IF(read_stat < 0) THEN
        RETURN 
      ENDIF  
        
      IF (trnspse == "T") THEN  
        DO row = 1,nrow
          READ(file_unit,*) (solution(map(col),row), col = 1,ncol)   
        ENDDO           
      ELSE
        DO row = 1,nrow
          READ(file_unit,*) (solution(map(row),col), col = 1,ncol)   
        ENDDO        
      ENDIF
     
      RETURN
      END SUBROUTINE read_solution_snap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE read_solution_full(out_direc,file_name,trnspse,t,solution,nsnap_read,lel2gel,last_snap)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc      
      CHARACTER(*), INTENT(IN) :: file_name
      CHARACTER(1), INTENT(IN) :: trnspse      
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: solution
      INTEGER, INTENT(INOUT) :: nsnap_read
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: lel2gel
      CHARACTER(1), INTENT(IN), OPTIONAL :: last_snap
      
      
      INTEGER :: el,n1,n2
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: nrow,ncol,nsnap
      INTEGER :: file_unit
      REAL(rp) :: ttemp
      INTEGER, DIMENSION(:), ALLOCATABLE :: map               
      CHARACTER(1) :: last_only
      
      CALL read_solution_header(out_direc,file_name,nrow,ncol,nsnap,file_unit)      
      
      IF (trnspse == "T") THEN
        n1 = ncol     
        n2 = nrow
      ELSE
        n1 = nrow
        n2 = ncol
      ENDIF
      
      IF (PRESENT(last_snap)) THEN
        last_only = last_snap
      ELSE
        last_only = "F"
      ENDIF
      
      
      IF (.not. ALLOCATED(solution)) THEN
        IF (last_only == "F") THEN
          ALLOCATE(solution(n1,n2,nsnap))
        ELSE
          ALLOCATE(solution(n1,n2,1))
        ENDIF
      ENDIF      
      
      IF (.not. ALLOCATED(t)) THEN
        ALLOCATE(t(nsnap))
      ENDIF
      
      ALLOCATE(map(n1))
      IF (PRESENT(lel2gel)) THEN
        DO el = 1,n1
          map(el) = lel2gel(el)
        ENDDO
      ELSE
        DO el = 1,n1
          map(el) = el
        ENDDO
      ENDIF   
      
      
      
      tstep = 0
      DO         
      
        tstep = tstep + 1   

        
        IF (last_only == "F") THEN
          CALL read_solution_snap(file_unit,nrow,ncol,trnspse,map,read_stat,ttemp,solution(:,:,tstep))
        ELSE
          CALL read_solution_snap(file_unit,nrow,ncol,trnspse,map,read_stat,ttemp,solution(:,:,1))
        ENDIF
        
        IF(read_stat < 0) THEN  
          tstep = tstep - 1
          EXIT 
        ENDIF
        
        t(tstep) = ttemp
        
        IF (tstep == nsnap_read) THEN        
          EXIT
        ENDIF
                
      ENDDO                   
      CLOSE(file_unit)
      
      IF (tstep < nsnap) THEN
        PRINT "(A)", "Number of solution snaps less than expected value"
      ENDIF
      
      nsnap_read = tstep

      RETURN
      END SUBROUTINE read_solution_full
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE copy_version_header(file_unit,filename)
      
      IMPLICIT NONE
            
      INTEGER, INTENT(IN) :: file_unit  
      CHARACTER(*), INTENT(IN) :: filename
      
      CHARACTER(200) :: line
                           
      OPEN(UNIT=9,file='PE0000/'//TRIM(ADJUSTL(filename)))
      
      DO
      
        READ(9,"(A)") line      
        
        IF (line(1:10) == "----------") THEN
          EXIT 
        ELSE 
          WRITE(file_unit,"(A)") line  
        ENDIF
      
      ENDDO 
   
      CLOSE(9)
      
      RETURN
      END SUBROUTINE copy_version_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE read_fort63(out_direc,t,eta,nsnap_read)
      
      IMPLICIT NONE            
      
      CHARACTER(*), INTENT(IN) :: out_direc
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: eta      
      INTEGER, INTENT(INOUT) :: nsnap_read
      
      INTEGER :: read_stat
      INTEGER :: k,nd,it,tstep
      INTEGER :: nsnap,nn,nskip,var_type
      REAL(rp) :: tf
      REAL(rp) :: ttemp
      CHARACTER(200) :: line      
      
      OPEN(UNIT=63,file=TRIM(ADJUSTL(out_direc)) // 'fort.63')      
      
      READ(63,*) line
      READ(63,*) nsnap,nn,tf,nskip,var_type
      
      
      IF (.not. ALLOCATED(eta)) THEN
        ALLOCATE(eta(nn,nsnap))
      ENDIF      
      
      IF (.not. ALLOCATED(t)) THEN
        ALLOCATE(t(nsnap))
      ENDIF      
      
      
      tstep = 0
      DO             
      
        READ(63,*,IOSTAT=read_stat) ttemp,it 
              
        IF(read_stat < 0) THEN  
          EXIT 
        ENDIF
        
        tstep = tstep + 1        
        t(tstep) = ttemp
                   
              
        DO nd = 1,nn        
          READ(63,*) k,eta(k,tstep)        
        ENDDO
        
        IF (tstep == nsnap_read) THEN        
          EXIT
        ENDIF           
            
      ENDDO
      
      CLOSE(63)
      
      IF (tstep < nsnap) THEN
        PRINT "(A)", "Number of solution snaps less than expected value"
      ENDIF
      
      nsnap_read = tstep      
      
      RETURN
      END SUBROUTINE read_fort63

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE read_fort64(out_direc,t,u,v,nsnap_read)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: u
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: v         
      INTEGER, INTENT(INOUT) :: nsnap_read      
      
      INTEGER :: read_stat
      INTEGER :: k,nd,it,tstep
      INTEGER :: nsnap,nn,nskip,var_type
      REAL(rp) :: tf
      REAL(rp) :: ttemp
      CHARACTER(200) :: line    
      
      OPEN(UNIT=64,file=TRIM(ADJUSTL(out_direc)) // 'fort.64')      
      
      READ(64,*) line
      READ(64,*) nsnap,nn,tf,nskip,var_type      
      
      
      
      IF (.not. ALLOCATED(u)) THEN
        ALLOCATE(u(nn,nsnap))
      ENDIF      
      
      IF (.not. ALLOCATED(v)) THEN
        ALLOCATE(v(nn,nsnap))
      ENDIF        
      
      IF (.not. ALLOCATED(t)) THEN
        ALLOCATE(t(nsnap))
      ENDIF      
      
      
      tstep = 0
      DO             
      
        READ(64,*,IOSTAT=read_stat) ttemp,it 
              
        IF(read_stat < 0) THEN  
          EXIT 
        ENDIF
        
        tstep = tstep + 1        
        t(tstep) = ttemp
                   
              
        DO nd = 1,nn        
          READ(64,*) k,u(k,tstep),v(k,tstep)
        ENDDO
        
        IF (tstep == nsnap_read) THEN        
          EXIT
        ENDIF           
            
      ENDDO      
      

      CLOSE(64)
      
      IF (tstep < nsnap) THEN
        PRINT "(A)", "Number of solution snaps less than expected value"
      ENDIF
      
      nsnap_read = tstep            
      RETURN
      END SUBROUTINE read_fort64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE read_write_output