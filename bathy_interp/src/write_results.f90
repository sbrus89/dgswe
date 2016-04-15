      SUBROUTINE write_results(mesh)

      USE globals, ONLY: grid

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      
      INTEGER :: sind,eind
      INTEGER :: et,nnd,el,i,j
      INTEGER :: flag
      CHARACTER(100) :: name
      CHARACTER(1) :: hbp


      sind = INDEX(ADJUSTL(TRIM(mesh%bathy_file)),"/",.true.)
      eind = INDEX(ADJUSTL(TRIM(mesh%grid_file)),".",.false.)
      
      
      name = ADJUSTL(TRIM(mesh%grid_file(1:eind-1))) 
      WRITE(hbp,"(I1)") mesh%hbp      
      
      OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(name)) // "_hbp" // hbp // "_interp.hb")            
      WRITE(13,"(2(I7,1x))") mesh%ne,mesh%hbp
      
      OPEN(UNIT = 14, FILE = TRIM(mesh%out_direc) // "elem_nodes.d")      
      WRITE(14,"(2(I7,1x))") mesh%ne,mesh%hbp      
      
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(6)          
        ENDIF    

        WRITE(13,"(2(I7),1x,60(e24.17,1x))") el,nnd,(mesh%elhb(i,el), i = 1,nnd)        
        WRITE(14,"(2(I7),1x,60(e24.17,1x))") el,nnd,(mesh%elhb(i,el), i = 1,nnd)
      ENDDO
      
      CLOSE(13)
      CLOSE(14)


       
      OPEN(unit=9, file=TRIM(mesh%out_direc) // "interp_nodes.d")
      WRITE(9,*) mesh%npts
      DO i = 1,mesh%npts
        WRITE(9,*) (mesh%hbxy(j,i), j = 1,3)
      ENDDO      
      CLOSE(9)
      
      OPEN(unit=10, file=TRIM(mesh%out_direc) // "boundary_nodes.d")
      WRITE(10,*) mesh%npts
      
      DO i = 1,mesh%nn
        WRITE(10,*) mesh%bnd_flag(i)
      ENDDO      
      
      DO i = 1,mesh%ned
        IF (mesh%ed_type(i) /= 0) THEN
          flag = 1
        ELSE 
          flag = 0
        ENDIF
        
        DO j = 1,mesh%hbp-1
          WRITE(10,*) flag
        ENDDO
      ENDDO      
      
      DO i = 1,mesh%ne
        WRITE(10,*) 0
      ENDDO
      CLOSE(10) 
      
      
      
      
      
      
      
      
      
      OPEN(unit=242,file=TRIM(mesh%out_direc) // 'bathy.d')
      DO el = 1,mesh%ne
      
        et = mesh%el_type(el)     
        IF (mod(et,2) == 1) THEN
          nnd = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(6)
        ENDIF
        
        DO i = 1,nnd
          WRITE(242,"(3(e24.17,1x))") mesh%elhbxy(i,el,1),mesh%elhbxy(i,el,2),mesh%elhb(i,el)
        ENDDO
      ENDDO
      CLOSE(242)      
      
      RETURN
      END SUBROUTINE write_results
