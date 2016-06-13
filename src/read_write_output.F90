      MODULE read_write_output
      
      USE globals, ONLY:rp      

      INTEGER :: unit_counter = 10
      
      CONTAINS            
      
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

      SUBROUTINE write_solution_snap(file_unit,ndof,ne,t_sol,solution,trnspse)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ndof      
      INTEGER, INTENT(IN) :: ne
      REAL(rp), INTENT(IN) :: t_sol
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: solution
      CHARACTER(1), INTENT(IN), OPTIONAL :: trnspse
      
      INTEGER :: dof
      INTEGER :: el
      INTEGER :: flip
      
      IF (PRESENT(trnspse)) THEN
        IF (trnspse == "T") THEN
          flip = 1
        ELSE
          flip = 0
        ENDIF
      ELSE
        flip = 0
      ENDIF
      
      WRITE(file_unit,"(e24.17)") t_sol
      IF (flip == 0) THEN
        DO dof = 1,ndof
          WRITE(file_unit,"(*(e24.17,1x))") (solution(el,dof), el=1,ne)
        ENDDO
      ELSE
        DO dof = 1,ndof
          WRITE(file_unit,"(*(e24.17,1x))") (solution(dof,el), el=1,ne)
        ENDDO      
      ENDIF
      
      
      RETURN
      END SUBROUTINE write_solution_snap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_stations_snap(file_unit,nsta,t_sta,stations)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: nsta
      REAL(rp), INTENT(IN) :: t_sta
      REAL(rp), DIMENSION(:), INTENT(IN) :: stations

      INTEGER :: sta
      
      WRITE(file_unit,"(e24.17)") t_sta
      DO sta = 1,nsta
        WRITE(file_unit,"(2x,e24.17)") stations(sta)
      ENDDO

      
      RETURN
      END SUBROUTINE write_stations_snap
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE write_solution_full(out_direc,file_name,ndof,ne,nsnap,t_sol,solution,post)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: ndof
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nsnap
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
        CALL write_solution_snap(file_unit,ndof,ne,t_sol(snap),solution(:,:,snap))
      ENDDO
      
      CLOSE(file_unit)
      
      RETURN
      END SUBROUTINE write_solution_full
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE write_stations_full(out_direc,file_name,nsta,nsnap,t_sta,stations,post)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(IN) :: nsta
      INTEGER, INTENT(IN) :: nsnap
      REAL(rp), DIMENSION(:), INTENT(IN) :: t_sta
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: stations
      LOGICAL, INTENT(IN), OPTIONAL :: post
      
      INTEGER :: file_unit
      INTEGER :: snap
      INTEGER :: ncol
      

      ncol = 1
      IF (PRESENT(post)) THEN ! if post, copy version information from PE0000
        CALL file_init(out_direc,file_name,nsta,ncol,nsnap,file_unit,post)
      ELSE
        CALL file_init(out_direc,file_name,nsta,ncol,nsnap,file_unit)      
      ENDIF
      
      DO snap = 1,nsnap
        CALL write_stations_snap(file_unit,nsta,t_sta(snap),stations(:,snap))
      ENDDO
            
      CLOSE(file_unit)
      
      RETURN
      END SUBROUTINE write_stations_full
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE read_header(out_direc,filename,nrow,ncol,nsnap,file_unit)
      
      IMPLICIT NONE
            
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(OUT) :: nrow
      INTEGER, INTENT(OUT) :: ncol
      INTEGER, INTENT(OUT) :: nsnap    
      INTEGER, INTENT(OUT) :: file_unit  
      
      CHARACTER(10) :: flag
                           
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
      END SUBROUTINE read_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE read_solution(out_direc,file_name,t,solution,lel2gel)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc      
      CHARACTER(*), INTENT(IN) :: file_name
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: solution
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: lel2gel
      
      INTEGER :: el,dof
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: globalize
      INTEGER :: ndof,ne,nsnap
      INTEGER :: file_unit
      REAL(rp) :: ttemp
      
      
      globalize = 0
      IF (PRESENT(lel2gel)) THEN
        globalize = 1
      ENDIF
      
      
      CALL read_header(out_direc,file_name,ndof,ne,nsnap,file_unit)          
      
      tstep = 0
      DO         
      
        tstep = tstep + 1   
        
        READ(file_unit,*,IOSTAT=read_stat) ttemp
        IF(read_stat < 0) THEN
          EXIT 
        ELSE  
          t(tstep) = ttemp
        ENDIF  
        
        IF (globalize) THEN
          DO dof = 1,ndof
            READ(file_unit,*) (solution(lel2gel(el),dof,tstep), el = 1,ne)   
          ENDDO           
        ELSE
          DO dof = 1,ndof
            READ(file_unit,*) (solution(el,dof,tstep), el = 1,ne)   
          ENDDO          
        ENDIF
                
      ENDDO                   
      CLOSE(file_unit)
      
      IF (tstep < nsnap) THEN
        PRINT "(A)", "Number of solution snaps less than expected value"
      ENDIF


      RETURN
      END SUBROUTINE read_solution
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      SUBROUTINE read_stations(out_direc,file_name,t,stations,sta_l2g)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: out_direc
      CHARACTER(*), INTENT(IN) :: file_name
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: t
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: stations
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: sta_l2g      
      
      INTEGER :: sta
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: globalize 
      INTEGER :: nsta,ncol,nsnap
      INTEGER :: file_unit
      REAL(rp) :: ttemp
      
      globalize = 0
      IF (PRESENT(sta_l2g)) THEN
        globalize = 1
      ENDIF      

      CALL read_header(out_direc,file_name,nsta,ncol,nsnap,file_unit) 
      
      tstep = 0
      DO         
      
        tstep = tstep + 1      
        
        READ(file_unit,*,IOSTAT=read_stat) ttemp
        IF(read_stat < 0) THEN
          EXIT 
        ELSE
          t(tstep) = ttemp
        ENDIF          
  
        IF (globalize) THEN
          DO sta = 1,nsta
            READ(file_unit,*) stations(sta_l2g(sta),tstep)   
          ENDDO
        ELSE  
          DO sta = 1,nsta
            READ(file_unit,*) stations(sta,tstep)   
          ENDDO          
        ENDIF
              
      ENDDO 
      CLOSE(file_unit)
      
      IF (tstep < nsnap) THEN
        PRINT "(A)", "Number of station snaps less than expected value"
      ENDIF      


      RETURN
      END SUBROUTINE read_stations
      
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

      END MODULE read_write_output