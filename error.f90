      PROGRAM error

      USE globals
      USE read_grid
      USE kdtree2_module
      USE evaluate
      USE basis

      IMPLICIT NONE

      INTEGER :: el,nd,pt,elc,elf,i,line,dof,srch
      INTEGER :: n1,n2,nvert,eln,clnd,etf,etc,n
      INTEGER :: ne,mndof
      INTEGER :: el_found
      INTEGER :: srchdp     
      REAL(pres) :: x(4),y(4),area,sarea,tol,found 
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: r,s,hb      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: xf,yf      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hc,Qxc,Qyc
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hf,Qxf,Qyf            
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: phi
      REAL(pres) :: HL2,QxL2,QyL2 
      REAL(pres) :: t
      TYPE(kdtree2), POINTER :: tree
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults
      
      tol = 1d-5
      srchdp = 4
            

      
      CALL read_input()
      CALL read_grids()
      
      i = 0
      DO el = 1,coarse%ne
        IF (coarse%bndel(el) == 1) THEN
          i = i + 1
          PRINT*, i,el
        ENDIF
      ENDDO
      
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
      CALL function_eval() ! evaluate basis and shape functions for fine grid quadrature points
      CALL detJ_eval()
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      ALLOCATE(xf(mnqpta),yf(mnqpta))
      ALLOCATE(r(mnqpta),s(mnqpta),hb(mnqpta))
      ALLOCATE(elf2elc(fine%ne))
      el_found = 0 
      
      ! Find coarse element each fine element is located in
      DO elf = 1,fine%ne
      
        etf = fine%el_type(elf)
        
        ! Compute the (x,y) coordinates of the first quadrature point (from fine element)
        xf(1) = 0d0
        yf(1) = 0d0
        DO nd = 1,fine%nnds(etf)
          xf(1) = xf(1) + fine%l(nd,1,etf)*fine%elxy(nd,elf,1)
          yf(1) = yf(1) + fine%l(nd,1,etf)*fine%elxy(nd,elf,2)
        ENDDO
      
        ! Find coarse node closest to quadrature point
        CALL kdtree2_n_nearest(tp=tree,qv=(/xf(1),yf(1)/),nn=srchdp,results=kdresults)
        
search: DO srch = 1,srchdp
          clnd = coarse%vxyn(kdresults(srch)%idx)      
        
          PRINT("(A,I5,A,I5,A,F11.4)"), "fine element #: ",elf, ",   closest coarse node: ", clnd, ",   distance: ", kdresults(1)%dis      
        
          ! Test elements associated with nearest node to see which coarse element contains the quadrature point
          found = 0  
    elem: DO elc = 1,coarse%nepn(clnd)
            eln = coarse%epn(elc,clnd)
          
            etc = coarse%el_type(eln)
            nvert = coarse%nverts(etc)                
          
            ! Compute the local (r,s) coarse element coordinates of the (x,y) quadrature point
            CALL newton(xf(1),yf(1),1,eln,r,s,hb)
          
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
          
            ! The quadrature point is in the element if the reference element area and sum of sub triangle are the same
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
            PRINT*, "element not found: going to next node"    ! Try elements around next closest node     
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
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      mndof = coarse%mndof
      ne = coarse%ne
      
      ALLOCATE(coarse%H(ne,mndof),coarse%Qx(ne,mndof),coarse%Qy(ne,mndof))
      ALLOCATE(coarse%phi(mndof,mnqpta,1))
      ALLOCATE(phi(mndof*mnqpta))
 
      mndof = fine%mndof
      ne = fine%ne
      
      ALLOCATE(fine%H(ne,mndof),fine%Qx(ne,mndof),fine%Qy(ne,mndof)) 
      
      ALLOCATE(Hc(mnqpta),Qxc(mnqpta),Qyc(mnqpta))
      ALLOCATE(Hf(mnqpta),Qxf(mnqpta),Qyf(mnqpta))
              
      
      OPEN(UNIT=11, FILE=trim(coarse%out_direc) //"solution_H.d")
      OPEN(UNIT=12, FILE=trim(coarse%out_direc) //"solution_Qx.d")     
      OPEN(UNIT=13, FILE=trim(coarse%out_direc) //"solution_Qy.d")      
      
      READ(11,*) 
      READ(12,*)
      READ(13,*)
      
      OPEN(UNIT=21, FILE=trim(fine%out_direc) //"solution_H.d")
      OPEN(UNIT=22, FILE=trim(fine%out_direc) //"solution_Qx.d")     
      OPEN(UNIT=23, FILE=trim(fine%out_direc) //"solution_Qy.d")      
      
      READ(21,*) 
      READ(22,*)
      READ(23,*)
      

      
      DO line = 1,lines+1
        READ(11,*) t
        READ(12,*) t
        READ(13,*) t
        
        READ(21,*) t
        READ(22,*) t
        READ(23,*) t        
        
        DO dof = 1,coarse%mndof
          READ(11,*) (coarse%H(el,dof), el = 1,coarse%ne)
          READ(12,*) (coarse%Qx(el,dof), el = 1,coarse%ne)
          READ(13,*) (coarse%Qy(el,dof), el = 1,coarse%ne)
        ENDDO    
        
        DO dof = 1,fine%mndof
          READ(21,*) (fine%H(el,dof), el = 1,fine%ne)
          READ(22,*) (fine%Qx(el,dof), el = 1,fine%ne)
          READ(23,*) (fine%Qy(el,dof), el = 1,fine%ne)
        ENDDO         
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)  
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      HL2 = 0d0
      QxL2 = 0d0
      QyL2 = 0d0
      
elemf:DO elf = 1,fine%ne        
        etf = fine%el_type(elf)
                 
        elc = elf2elc(elf)    
        
        IF (coarse%bndel(elc) == 1) THEN
          CYCLE elemf
        ENDIF
        
        etc = coarse%el_type(elc)  
        
        xf = 0d0
        yf = 0d0
        DO pt = 1,nqpta(etf)
          DO nd = 1,fine%nnds(etf)
            xf(pt) = xf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,1)
            yf(pt) = yf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,2)
          ENDDO        
        ENDDO
        
        CALL newton(xf,yf,nqpta(etf),elc,r,s,hb)
                 
        IF (mod(etc,2) == 1) THEN
          CALL tri_basis(coarse%p,coarse%ndof(etc),nqpta(etf),r,s,phi)       
        ELSE IF (mod(etc,2) == 0) THEN
          CALL quad_basis(coarse%p,coarse%ndof(etc),nqpta(etf),r,s,phi)
        ENDIF          
          
        DO pt = 1,nqpta(etf)
          DO dof = 1,coarse%ndof(etc)
             i = (dof-1)*nqpta(etf) + pt
            coarse%phi(dof,pt,1) = phi(i)
          ENDDO
        ENDDO
          
        Hc = 0d0
        Qxc = 0d0
        Qyc = 0d0
        DO pt = 1,nqpta(etf)
          DO dof = 1,coarse%ndof(etc)
            Hc(pt) = Hc(pt) + coarse%H(elc,dof)*coarse%phi(dof,pt,1)
            Qxc(pt) = Qxc(pt) + coarse%Qx(elc,dof)*coarse%phi(dof,pt,1)            
            Qyc(pt) = Qyc(pt) + coarse%Qy(elc,dof)*coarse%phi(dof,pt,1)
          ENDDO                
        ENDDO
        
        Hf = 0d0
        Qxf = 0d0
        Qyf = 0d0
        DO pt = 1,nqpta(etf)
          DO dof = 1,fine%ndof(etf)
            Hf(pt) = Hf(pt) + fine%H(elf,dof)*fine%phi(dof,pt,etf)
            Qxf(pt) = Qxf(pt) + fine%Qx(elf,dof)*fine%phi(dof,pt,etf)            
            Qyf(pt) = Qyf(pt) + fine%Qy(elf,dof)*fine%phi(dof,pt,etf)
          ENDDO                
        ENDDO       
        
        DO pt = 1,nqpta(etf)
          HL2 = HL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Hc(pt)-Hf(pt))**2
          QxL2 = QxL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Qxc(pt)-Qxf(pt))**2
          QyL2 = QyL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Qyc(pt)-Qyf(pt))**2
        ENDDO
        
      ENDDO elemf    
      
      HL2 = sqrt(HL2)
      QxL2 = sqrt(QxL2)
      QyL2 = sqrt(QyL2)
      
      PRINT*, HL2, QxL2,QyL2

          
      END PROGRAM error