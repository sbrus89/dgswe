      PROGRAM error

      USE globals
      USE read_grid
      USE kdtree2_module
      USE evaluate
      USE basis

      IMPLICIT NONE

      INTEGER :: el,nd,pt,elc,elf,i,line,dof,srch
      INTEGER :: n1,n2,nvert,eln,clnd,etf,etc,n
      INTEGER :: el_found
      INTEGER :: srchdp
      REAL(pres) :: rsre(2,4,4),x(4),y(4),area,sarea,tol,found 
      REAL(pres) :: xyf(2)
      REAL(pres) :: r(1),s(1),hb
      REAL(pres) :: Hpt,Qxpt,Qypt
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: phi
      REAL(pres) :: t
      TYPE(kdtree2), POINTER :: tree
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults
      
      tol = 1d-5
      srchdp = 4
            

      
      CALL read_input()
      CALL read_grids()
      
                   
      ! Build kd-tree
      IF (coarse%ne < srchdp) THEN
        srchdp = coarse%ne
      ENDIF      
      ALLOCATE(kdresults(srchdp))            
      tree => kdtree2_create(coarse%vxy, rearrange=.true., sort=.true.)
      
      
      CALL vandermonde(coarse)
      CALL vandermonde(fine)
            
      CALL re_vert(coarse)
      CALL re_vert(fine)
      
      CALL area_qpts()      
      CALL l_eval()
      
      
      ALLOCATE(elf2elc(fine%ne))
      el_found = 0 
      
      ! Find element each station is located in and determine local coordinates of station
      DO elf = 1,fine%ne
      
        etf = fine%el_type(elf)
      
        xyf = 0d0
        DO nd = 1,fine%nnds(etf)
          xyf(1) = xyf(1) + fine%l(nd,1,etf)*fine%elxy(nd,elf,1)
          xyf(2) = xyf(2) + fine%l(nd,1,etf)*fine%elxy(nd,elf,2)
        ENDDO
      
        ! Find node closest to station
        CALL kdtree2_n_nearest(tp=tree,qv=xyf,nn=srchdp,results=kdresults)
        
search: DO srch = 1,srchdp
          clnd = coarse%vxyn(kdresults(srch)%idx)      
        
          PRINT("(A,I5,A,I5,A,F11.4)"), "fine element #: ",elf, ",   closest coarse node: ", clnd, ",   distance: ", kdresults(1)%dis      
        
          ! Test elements associated with nearest node to see which element station is located in
          found = 0  
    elem: DO elc = 1,coarse%nepn(clnd)
            eln = coarse%epn(elc,clnd)
          
            etc = coarse%el_type(eln)
            nvert = coarse%nverts(etc)                
          
            ! Compute the local (r,s) coordinates of the (x,y) station location
            CALL newton(xyf(1),xyf(2),eln,r,s,hb)
          
            ! Find reference element area
            IF (mod(etc,2) == 1) THEN
              area = 2d0
            ELSE IF (mod(etc,2) == 0) THEN
              area = 4d0
            ENDIF          
          
            ! Compute sum of sub-triangle areas
            sarea = 0d0
            DO i = 1,nvert
              n1 = mod(i+0,nvert)+1
              n2 = mod(i+1,nvert)+1           

              x(1) = coarse%rsre(1,n1,etc)
              y(1) = coarse%rsre(2,n1,etc)
            
              x(2) = coarse%rsre(1,n2,etc)
              y(2) = coarse% rsre(2,n2,etc)
            
              x(3) = r(1)
              y(3) = s(1)
            
              sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
            ENDDO
          
            PRINT("(A,I5,A,F20.15,A,F20.15)"), "   testing: ", eln, "   area = ",area, "   sarea = ", sarea
            PRINT*, " "
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            IF (abs(area - sarea) < tol) THEN
              PRINT("(A,I5)"), '   element found ', eln            
!               PRINT("(A,I5)"), achar(27)//'[95m    element found '//achar(27)//'[0m', eln
            
              elf2elc(elf) = eln            
            
              found = 1
              el_found = el_found + 1
                        
            
              EXIT elem                        
            ENDIF
          
          ENDDO elem
        
          IF (found == 0) THEN
            PRINT*, "element not found: going to next node"         
          ELSE   
            EXIT search
          ENDIF
        
          PRINT*, " " 
        
        ENDDO search
        
        IF (found == 0) THEN
          PRINT*,  achar(27)//"[31m ERROR: ELEMENT NOT FOUND"//achar(27)//'[0m'
        ENDIF        
      
      ENDDO
      
      PRINT*, "Missing elements = ", fine%ne-el_found
      PRINT*, " "
      
      IF (fine%ne-el_found /= 0) THEN
        STOP
      ENDIF
!         
! 
!       
!       ALLOCATE(H(ne,mndof),Qx(ne,mndof),Qy(ne,mndof))
!       ALLOCATE(phi(mndof))
!               
!       
!       OPEN(UNIT=63, FILE=trim(out_direc) //"solution_H.d")
!       OPEN(UNIT=641,FILE=trim(out_direc) //"solution_Qx.d")     
!       OPEN(UNIT=642,FILE=trim(out_direc) //"solution_Qy.d")      
!       
!       READ(63,*) 
!       READ(641,*)
!       READ(642,*)
! 
!       
!       DO line = 1,lines+1
!         READ(63,*) t
!         READ(641,*) t
!         READ(642,*) t
!         
!         PRINT("(A,I5)"), "Writing time snap: ", line
!         
!         DO dof = 1,mndof
!           READ(63,*) (H(el,dof), el = 1,ne)
!           READ(641,*) (Qx(el,dof), el = 1,ne)
!           READ(642,*) (Qy(el,dof), el = 1,ne)
!         ENDDO
!                
!                 
!         DO sta = 1,nsta
!                  
!           eln = elsta(sta)                 
!           et = el_type(eln)                 
!                  
!           IF (mod(et,2) == 1) THEN
!             CALL tri_basis(p,ndof(et),1,rssta(sta,1),rssta(sta,2),phi)       
!           ELSE IF (mod(et,2) == 0) THEN
!             CALL quad_basis(p,ndof(et),1,rssta(sta,1),rssta(sta,2),phi)
!           ENDIF          
!           
!           Hpt = 0d0
!           Qxpt = 0d0
!           Qypt = 0d0
!           DO dof = 1,ndof(et)
!             Hpt = Hpt + H(eln,dof)*phi(dof)
!             Qxpt = Qxpt + Qx(eln,dof)*phi(dof)            
!             Qypt = Qypt + Qy(eln,dof)*phi(dof)
!           ENDDO                
!           
!         ENDDO
!         
!       ENDDO
!       
!       CLOSE(63)
!       CLOSE(641)
!       CLOSE(642)
          
      END PROGRAM error