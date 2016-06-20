      MODULE plot_mod
      
      USE globals, ONLY: rp
      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      
      SUBROUTINE write_psheader()
      
      OPEN(UNIT=100,FILE='out.ps')      
      WRITE(100,"(A)") "%!PS-Adobe-3.0"
      
      WRITE(100,"(A)") "/draw-element {"
      WRITE(100,"(A)") "moveto"
      WRITE(100,"(A)") "lineto"
      WRITE(100,"(A)") "lineto"
      WRITE(100,"(A)") "closepath"
      WRITE(100,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(100,"(A)") "stroke"      
      WRITE(100,"(A)") "} def" 
      
      WRITE(100,"(A)") "/trifill {"
      WRITE(100,"(A)") "/A exch def"
      WRITE(100,"(A)") "/B exch def"
      WRITE(100,"(A)") "/C exch def"     
      WRITE(100,"(A)") "/ds ["
      WRITE(100,"(A)") "  0 A aload pop" 
      WRITE(100,"(A)") "  0 B aload pop"
      WRITE(100,"(A)") "  0 C aload pop"
      WRITE(100,"(A)") "] def"
      WRITE(100,"(A)") "newpath"
      WRITE(100,"(A)") "<<"
      WRITE(100,"(A)") "  /ShadingType 4"
      WRITE(100,"(A)") "  /ColorSpace [ /DeviceRGB ]"
      WRITE(100,"(A)") "  /DataSource ds" 
!       WRITE(100,"(A)")   " /AntiAlias true"
      WRITE(100,"(A)") ">>"
      WRITE(100,"(A)") "shfill"
      WRITE(100,"(A)") "} def"        
      
      END SUBROUTINE write_psheader
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      
      SUBROUTINE close_ps()

      WRITE(100,"(A)") "showpage"
      CLOSE(100)            
      
      END SUBROUTINE close_ps     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE plot_contours(nplt,ntri,rect,ne,el_type,el_in,xy,sol_val)
      
      IMPLICIT NONE
      
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      INTEGER, DIMENSION(:), INTENT(IN) :: ntri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xy  
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sol_val
      
      INTEGER :: i,j,v
      INTEGER :: el,nd,dof,lev,tri
      INTEGER :: et,nv
      INTEGER :: nlev      
      REAL(rp) :: sol_min,sol_max,sol_lev
      REAL(rp) :: dc
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors     
      REAL(rp) :: color_val(3)
      REAL(rp) :: cmax



      sol_min = 999d0
      sol_max = -999d0
      
      DO el = 1,ne
        et = el_type(el)
        nv = nplt(et)
        DO nd = 1,nv
          IF (sol_val(nd,el) < sol_min) THEN
            sol_min = sol_val(nd,el)
          ENDIF
        
          IF (sol_val(nd,el) > sol_max) THEN
            sol_max = sol_val(nd,el)
          ENDIF          
        ENDDO
        
      ENDDO
      
      PRINT*, sol_min,sol_max
      
  

      OPEN(UNIT=101,FILE="default2.cmap")
      READ(101,*) nlev
      ALLOCATE(colors(nlev,3))
      DO lev = 1,nlev
        READ(101,*) (colors(lev,j), j=1,3)
      ENDDO
      CLOSE(101)
      
      cmax = MAXVAL(colors)
      IF (cmax > 1d0+1d-12) THEN
        colors = colors/255d0
      ENDIF
      
      dc = (sol_max-sol_min)/real(nlev-1,rp)
            
      
 elem:DO el = 1,ne
        et = el_type(el)
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        DO tri = 1,ntri(et)

        DO v = 1,3  
   
          nd = rect(v,tri,et)
          
          sol_lev = sol_min
          lev = nlev
  levels: DO i = 1,nlev-1
            IF ((sol_val(nd,el) >= sol_lev) .and. (sol_val(nd,el) < sol_lev+dc)) THEN
              lev = i
              CALL interp_colors(lev,sol_lev,dc,colors,sol_val(nd,el),color_val)
              EXIT levels
            ENDIF
            sol_lev = sol_lev + dc
          ENDDO levels
          
          IF (sol_val(nd,el) <= sol_min) THEN
           lev = 1
         ELSE IF (sol_val(nd,el) > sol_max) THEN
           lev = nlev
         ENDIF          

          WRITE(100,"(A,5(F9.5,1x),A)") "[",xy(nd,el,1),xy(nd,el,2),color_val(1),color_val(2),color_val(3),"]"  
        ENDDO        

        WRITE(100,"(A)") "trifill"        
        
        ENDDO
      ENDDO elem


      END SUBROUTINE plot_contours
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_mesh(ne,nverts,el_type,el_in,xy,ect)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect

      INTEGER :: el,nd
      INTEGER :: et
      
 elem:DO el = 1,ne
        et = el_type(el)
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        DO nd = 1,nverts(et)
          WRITE(100,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
        ENDDO        

        WRITE(100,"(A)") "draw-element"
      ENDDO elem

      END SUBROUTINE plot_mesh
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE interp_colors(lev,sol_lev,dc,colors,sol_val,color_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: lev
      REAL(rp), INTENT(IN) :: sol_lev
      REAL(rp), INTENT(IN) :: dc
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: colors
      REAL(rp), INTENT(IN) :: sol_val
      REAL(rp), DIMENSION(:), INTENT(OUT) :: color_val
      
      INTEGER :: i
      REAL(rp) :: l1,l2
      
      l1 = (sol_val-(sol_lev+dc))/(-dc)
      l2 = (sol_val-sol_lev)/dc
      
      DO i = 1,3
        color_val(i) = l1*colors(lev,i) + l2*colors(lev+1,i)
      ENDDO      
      
      
!       DO i = 1,3
!         color_val(i) = colors(lev,i)
!       ENDDO
!       
      RETURN
      END SUBROUTINE interp_colors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      END MODULE plot_mod