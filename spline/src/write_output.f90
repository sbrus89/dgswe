
      
      SUBROUTINE write_grid(mesh)

      USE globals, ONLY: grid,ctp

      IMPLICIT NONE
          
      INTEGER :: nd,el,bnd
      TYPE(grid) :: mesh
      
      OPEN(UNIT=14,FILE='nodes.out')
      WRITE(14,*) "spline output nodes"
      WRITE(14,*) mesh%ne,mesh%nn
      DO nd = 1,mesh%nn
        WRITE(14,"(I7,3(e24.16))") nd,mesh%xy(1,nd),mesh%xy(2,nd),mesh%depth(nd)
      ENDDO
      DO el = 1,mesh%ne
        WRITE(14,"(20(I9))") el,mesh%nelnds(el), (mesh%ect(nd,el),nd=1,mesh%mnelnds)
      ENDDO
      
      WRITE(14,"(I9,A)") mesh%nope, "   = Number of open boundaries"
      WRITE(14,"(I9,A)") mesh%neta, "   = Total number of open boundary nodes"
      DO bnd = 1,mesh%nope
        WRITE(14,"(2(I9))") mesh%obseg(bnd), 0
        DO nd = 1,mesh%obseg(bnd)
          WRITE(14,"(I9)") mesh%obnds(nd,bnd)
        ENDDO
      ENDDO
      
      WRITE(14,"(I9,A)") mesh%nbou, "   = Number of land boundaries"
      WRITE(14,"(I9,A)") mesh%nvel, "   = Total number of land boundary nodes"
      DO bnd = 1,mesh%nbou
        WRITE(14,"(2(I9))") mesh%fbseg(1,bnd), mesh%fbseg(2,bnd)
        DO nd = 1,mesh%fbseg(1,bnd)
          WRITE(14,"(I9)") mesh%fbnds(nd,bnd)
        ENDDO
      ENDDO
      
      CLOSE(14)

      RETURN
      END SUBROUTINE write_grid
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE write_spline(n)
      
      USE globals, ONLY: rp,ax,bx,cx,dx,ay,by,cy,dy,dt
      
      IMPLICIT NONE
           
      INTEGER :: n      
      INTEGER :: i
      REAL(rp) :: t
      
      
      WRITE(30,*) n

      t = 0d0
      DO i = 1,n-1        
        WRITE(30,"(9(E25.12,1x))"), ax(i),bx(i),cx(i),dx(i),ay(i),by(i),cy(i),dy(i),t
        t = t + dt(i)
      ENDDO
      WRITE(30,"(9(E25.12,1x))"), ax(n),0d0,0d0,0d0,ay(n),0d0,0d0,0d0,t

      
      
      RETURN
      END SUBROUTINE write_spline
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      