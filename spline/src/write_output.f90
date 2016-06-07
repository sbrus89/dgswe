
      
      SUBROUTINE write_grid(mesh,base_name)

      USE globals, ONLY: grid,nverts

      IMPLICIT NONE
          
      INTEGER :: nd,el,bnd,et,nv,btype,nbseg
      TYPE(grid) :: mesh
      CHARACTER(100) :: base_name
      
      OPEN(UNIT=14,FILE='nodes.out')
      WRITE(14,"(A)") "spline output nodes based on "//base_name
      WRITE(14,*) mesh%ne,mesh%nn
      DO nd = 1,mesh%nn
        WRITE(14,"(I7,3(e24.16))") nd,mesh%xy(1,nd),mesh%xy(2,nd),mesh%depth(nd)
      ENDDO
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        nv = nverts(et)
        WRITE(14,"(20(I9))") el,nv, (mesh%ect(nd,el),nd=1,nv)
      ENDDO
      
      WRITE(14,"(I9,A)") mesh%nope, "   = Number of open boundaries"
      WRITE(14,"(I9,A)") mesh%neta, "   = Total number of open boundary nodes"
      DO bnd = 1,mesh%nope
        WRITE(14,"(I9)") mesh%obseg(bnd)
        DO nd = 1,mesh%obseg(bnd)
          WRITE(14,"(I9)") mesh%obnds(nd,bnd)
        ENDDO
      ENDDO
      
      WRITE(14,"(I9,A)") mesh%nbou, "   = Number of land boundaries"
      WRITE(14,"(I9,A)") mesh%nvel, "   = Total number of land boundary nodes"
      DO bnd = 1,mesh%nbou
        btype = mesh%fbseg(2,bnd)
        nbseg = mesh%fbseg(1,bnd)
        
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          nbseg = nbseg - 1
        ENDIF        
        
        WRITE(14,"(2(I9))") nbseg, btype
        DO nd = 1,nbseg
          WRITE(14,"(I9)") mesh%fbnds(nd,bnd)
        ENDDO        
      ENDDO
      
      CLOSE(14)

      RETURN
      END SUBROUTINE write_grid
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE write_spline(n,bou,ndnums)
      
      USE globals, ONLY: rp,ax,bx,cx,dx,ay,by,cy,dy,dt
      
      IMPLICIT NONE
           
      INTEGER :: n,bou      
      INTEGER :: i
      INTEGER, DIMENSION(n) :: ndnums
      REAL(rp) :: t
      
      
      WRITE(30,*) n

      t = 0d0
      DO i = 1,n-1        
        WRITE(30,"(9(E25.12,1x),I7)"), ax(i,bou),bx(i,bou),cx(i,bou),dx(i,bou),ay(i,bou),by(i,bou),cy(i,bou),dy(i,bou),t,ndnums(i)
        t = t + dt(i,bou)
      ENDDO
      WRITE(30,"(9(E25.12,1x),I7)"), ax(n,bou),0d0,0d0,0d0,ay(n,bou),0d0,0d0,0d0,t,ndnums(n)

      
      
      RETURN
      END SUBROUTINE write_spline
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      