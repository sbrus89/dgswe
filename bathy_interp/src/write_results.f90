      SUBROUTINE write_results(eval,base)

      USE globals, ONLY: grid
      USE version, ONLY: version_information

      IMPLICIT NONE
      
      TYPE(grid) :: eval,base
      
      INTEGER :: sind,eind
      INTEGER :: et,nnd,el,i,j
      INTEGER :: flag
      CHARACTER(100) :: name
      CHARACTER(1) :: hbp_char


!       sind = INDEX(ADJUSTL(TRIM(eval%bathy_file)),"/",.true.)
!       eind = INDEX(ADJUSTL(TRIM(eval%grid_file)),".",.false.)
!       
!       
!       name = ADJUSTL(TRIM(eval%grid_file(1:eind-1))) 
!       WRITE(hbp,"(I1)") eval%hbp      
!       
!       OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(name)) // "_hbp" // hbp // "_interp.hb") 
!       WRITE(13,"(A,I5)") " base grid: " // ADJUSTL(TRIM(base%grid_file)) // "    base hbp: ", base%hbp       
!       WRITE(13,"(2(I7,1x))") eval%ne,eval%hbp
!       
!       OPEN(UNIT = 14, FILE = TRIM(eval%out_direc) // "elem_nodes.d")      
!       WRITE(14,"(2(I7,1x))") eval%ne,eval%hbp      
!       
!       DO el = 1,eval%ne
!         et = eval%el_type(el)
!         
!         IF (mod(et,2) == 1) THEN   
!           nnd = eval%nnds(5)
!         ELSE IF (mod(et,2) == 0) THEN
!           nnd = eval%nnds(6)          
!         ENDIF    
! 
!         WRITE(13,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elhb(i,el), i = 1,nnd)        
!         WRITE(14,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elhb(i,el), i = 1,nnd)
!       ENDDO
!       
!       CLOSE(13)
!       CLOSE(14)



      eind = INDEX(ADJUSTL(TRIM(eval%grid_file)),".",.false.)
      name = ADJUSTL(TRIM(eval%grid_file(1:eind-1)))           
      WRITE(hbp_char,"(I1)") eval%hbp
      
      OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(name)) // "_hbp" // hbp_char // "_interp.hb") 
      
      CALL version_information(13)
      
      WRITE(13,"(A)") "-----------------------------------------------------------------------"           
      
      CALL write_input(13)           
      
      WRITE(13,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"         

      WRITE(13,"(2(I7,1x),A,I8)") eval%ne,eval%hbp      
      
      DO el = 1,eval%ne
        et = eval%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = eval%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = eval%nnds(6)          
        ENDIF    

        WRITE(13,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elhb(i,el), i = 1,nnd)        
      ENDDO      
      
      CLOSE(13)



      
      OPEN(UNIT = 14, FILE = TRIM(eval%out_direc) // "elem_nodes.d")      
      WRITE(14,"(2(I7,1x))") eval%ne,eval%hbp      
      
      DO el = 1,eval%ne
        et = eval%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = eval%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = eval%nnds(6)          
        ENDIF    

        WRITE(14,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elhb(i,el), i = 1,nnd)
      ENDDO
      
      CLOSE(14)

       
      OPEN(unit=9, file=TRIM(eval%out_direc) // "interp_nodes.d")
      WRITE(9,*) eval%npts
      DO i = 1,eval%npts
        WRITE(9,*) (eval%hbxy(j,i), j = 1,3)
      ENDDO      
      CLOSE(9)
      
      OPEN(unit=10, file=TRIM(eval%out_direc) // "boundary_nodes.d")
      WRITE(10,*) eval%npts
      
      DO i = 1,eval%nn
        WRITE(10,*) eval%bnd_flag(i)
      ENDDO      
      
      DO i = 1,eval%ned
        IF (eval%ed_type(i) /= 0) THEN
          flag = 1
        ELSE 
          flag = 0
        ENDIF
        
        DO j = 1,eval%hbp-1
          WRITE(10,*) flag
        ENDDO
      ENDDO      
      
      DO i = 1,eval%ne
        WRITE(10,*) 0
      ENDDO
      CLOSE(10) 
      
      
      
      
      
      
      
      
      
      OPEN(unit=242,file=TRIM(eval%out_direc) // 'bathy.d')
      DO el = 1,eval%ne
      
        et = eval%el_type(el)     
        IF (mod(et,2) == 1) THEN
          nnd = eval%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = eval%nnds(6)
        ENDIF
        
        DO i = 1,nnd
          WRITE(242,"(3(e24.17,1x))") eval%elhbxy(i,el,1),eval%elhbxy(i,el,2),eval%elhb(i,el)
        ENDDO
      ENDDO
      CLOSE(242)      
      
      RETURN
      END SUBROUTINE write_results
