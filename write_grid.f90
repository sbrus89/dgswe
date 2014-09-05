      SUBROUTINE write_grid()

      USE globals, ONLY: pres,ne,nn,ctp,xy,ect,depth,nelnds, &
                         nope,neta,obnds,obseg,nbou,nvel,fbseg,fbnds, &
                         mnelnds

      IMPLICIT NONE
      
      INTEGER :: nd,el,bnd
      
      OPEN(UNIT=14,FILE='nodes.out')
      WRITE(14,*) "spline output nodes"
      WRITE(14,*) ne,nn
      DO nd = 1,nn
        WRITE(14,"(I7,3(e24.16))") nd,xy(1,nd),xy(2,nd),depth(nd)
      ENDDO
      DO el = 1,ne
        WRITE(14,"(20(I9))") el,nelnds(el), (ect(nd,el),nd=1,mnelnds)
      ENDDO
      
      WRITE(14,"(I9,A)") nope, "   = Number of open boundaries"
      WRITE(14,"(I9,A)") neta, "   = Total number of open boundary nodes"
      DO bnd = 1,nope
        WRITE(14,"(2(I9))") obseg(bnd), 0
        DO nd = 1,obseg(bnd)
          WRITE(14,"(I9)") obnds(nd,bnd)
        ENDDO
      ENDDO
      
      WRITE(14,"(I9,A)") nbou, "   = Number of land boundaries"
      WRITE(14,"(I9,A)") nvel, "   = Total number of land boundary nodes"
      DO bnd = 1,nbou
        WRITE(14,"(2(I9))") fbseg(1,bnd), fbseg(2,bnd)
        DO nd = 1,fbseg(1,bnd)
          WRITE(14,"(I9)") fbnds(nd,bnd)
        ENDDO
      ENDDO
      
      CLOSE(14)

      RETURN
      END SUBROUTINE write_grid