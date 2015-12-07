      SUBROUTINE write_results(mesh)

      USE globals, ONLY: grid

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      
      INTEGER :: sind
      INTEGER :: et,nnd,el,i,j
      CHARACTER(100) :: out_direc


      sind = INDEX(ADJUSTL(TRIM(mesh%bathy_file)),"/",.true.)
      
      out_direc = ADJUSTL(TRIM(mesh%bathy_file(1:sind)))     
      
      OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(mesh%bathy_file)))            
      WRITE(13,"(2(I7,1x))") mesh%ne,mesh%hbp
      
      OPEN(UNIT = 14, FILE = TRIM(out_direc) // "elem_nodes.d")      
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


       
      OPEN(unit=9, file=TRIM(out_direc) // "interp_nodes.d")
      WRITE(9,*) mesh%npts
      DO i = 1,mesh%npts
        WRITE(9,*) (mesh%hbxy(j,i), j = 1,3)
      ENDDO      
      CLOSE(9)
      
      RETURN
      END SUBROUTINE write_results
