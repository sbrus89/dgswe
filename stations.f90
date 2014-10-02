      PROGRAM stations

      USE globals
      USE kdtree2_module
      USE evaluate
      USE basis

      IMPLICIT NONE

      INTEGER :: el,nd,pt,sta,i,line,dof,srch
      INTEGER :: n1,n2,nvert,eln,clnd,et,n
      INTEGER :: el_found
      INTEGER :: mnepn
      INTEGER :: srchdp
      REAL(pres) :: rsre(2,4,4),x(4),y(4),area,sarea,tol,found 
      REAL(pres) :: r(1),s(1),hb
      REAL(pres) :: Hsta,Qxsta,Qysta
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: phi
      REAL(pres) :: t
      TYPE(kdtree2), POINTER :: tree
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults
      
      tol = 1d-5
      srchdp = 4
            

      
      CALL read_input()
      CALL read_grid()
      CALL read_stations()
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)*(p+1)
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)
      mndof = maxval(ndof)
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp        
      
      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      mnnds = maxval(nnds)          
      
      
      ! Build kd-tree
      IF (ne < srchdp) THEN
        srchdp = ne
      ENDIF      
      ALLOCATE(kdresults(srchdp))            
      tree => kdtree2_create(vxy, rearrange=.true., sort=.true.)
      
      
      
      
      
      ! Find elements associated with each node            
      ALLOCATE(nepn(nn))
      
      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
!         n = nelnds(el)
        DO nd = 1,nvert
          n1 = vct(nd,el)
          nepn(n1) = nepn(n1) + 1
        ENDDO
      ENDDO
      
      mnepn = maxval(nepn)
      
      ALLOCATE(epn(mnepn,nn))
      
      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
!         n = nelnds(el)
        DO nd = 1,nvert
          n1 = vct(nd,el)
          nepn(n1) = nepn(n1) + 1
          epn(nepn(n1),n1) = el
        ENDDO
      ENDDO      
      
      
      
      
      CALL vandermonde()  
      
      DO et = 1,nel_type
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(1),nnds(1),rsre(1,:,et),rsre(2,:,et))
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(2),nnds(2),rsre(1,:,et),rsre(2,:,et))
        ENDIF
        
!         PRINT*, (rsre(1,pt,et), pt = 1,4)
!         PRINT*, (rsre(2,pt,et), pt = 1,4)
!         
!         PRINT*, " "
      ENDDO
      
      
      ALLOCATE(rssta(nsta,2))
      ALLOCATE(elsta(nsta),hbsta(nsta))
      el_found = 0 
      
      ! Find element each station is located in and determine local coordinates of station
      DO sta = 1,nsta
      
        ! Find node closest to station
        CALL kdtree2_n_nearest(tp=tree,qv=xysta(:,sta),nn=srchdp,results=kdresults)
        
search: DO srch = 1,srchdp
          clnd = vxyn(kdresults(srch)%idx)      
        
          PRINT("(A,I5,A,I5,A,F11.4)"), "station #: ",sta, ",   closest node: ", clnd, ",   distance: ", kdresults(1)%dis      
        
          ! Test elements associated with nearest node to see which element station is located in
          found = 0  
    elem: DO el = 1,nepn(clnd)
            eln = epn(el,clnd)
          
            et = el_type(eln)
            nvert = nverts(et)      
          
          
            ! Compute the local (r,s) coordinates of the (x,y) station location
            CALL newton(xysta(1,sta),xysta(2,sta),eln,r,s,hb)
          
            ! Find reference element area
            IF (mod(et,2) == 1) THEN
              area = 2d0
            ELSE IF (mod(et,2) == 0) THEN
              area = 4d0
            ENDIF          
          
            ! Compute sum of sub-triangle areas
            sarea = 0d0
            DO i = 1,nvert
              n1 = mod(i+0,nvert)+1
              n2 = mod(i+1,nvert)+1           

              x(1) = rsre(1,n1,et)
              y(1) = rsre(2,n1,et)
            
              x(2) = rsre(1,n2,et)
              y(2) = rsre(2,n2,et)
            
              x(3) = r(1)
              y(3) = s(1)
            
              sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
            ENDDO
          
            PRINT("(A,I5,A,F20.15,A,F20.15)"), "   testing: ", eln, "   area = ",area, "   sarea = ", sarea
            PRINT*, " "
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            IF (abs(area - sarea) < tol) THEN
              PRINT("(A,I5)"), "   element found", eln
            
              elsta(sta) = eln
            
              rssta(sta,1) = r(1)
              rssta(sta,2) = s(1)
              hbsta(sta) = hb
            
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
          PRINT*, "ERROR: ELEMENT NOT FOUND"
        ENDIF        
      
      ENDDO
      
      PRINT*, "Missing elements = ", nsta-el_found
      PRINT*, " "
      
      IF (nsta-el_found /= 0) THEN
        STOP
      ENDIF
        

      
      ALLOCATE(H(ne,mndof),Qx(ne,mndof),Qy(ne,mndof))
      ALLOCATE(phi(mndof))
      
      OPEN(UNIT=61,FILE=trim(out_direc) //"station_H.d")
      OPEN(UNIT=621,FILE=trim(out_direc) //"station_Qx.d")     
      OPEN(UNIT=622,FILE=trim(out_direc) //"station_Qy.d")           
      
      OPEN(UNIT=63, FILE=trim(out_direc) //"solution_H.d")
      OPEN(UNIT=641,FILE=trim(out_direc) //"solution_Qx.d")     
      OPEN(UNIT=642,FILE=trim(out_direc) //"solution_Qy.d")      
      
      READ(63,*) 
      READ(641,*)
      READ(642,*)
      
      WRITE(61,*) nsta
      WRITE(621,*) nsta
      WRITE(622,*) nsta
      
      DO line = 1,lines+1
        READ(63,*) t
        READ(641,*) t
        READ(642,*) t
        
        PRINT("(A,I5)"), "Writing time snap: ", line
        
        DO dof = 1,mndof
          READ(63,*) (H(el,dof), el = 1,ne)
          READ(641,*) (Qx(el,dof), el = 1,ne)
          READ(642,*) (Qy(el,dof), el = 1,ne)
        ENDDO
        
        WRITE(61,*) t
        WRITE(621,*) t
        WRITE(622,*) t           
                
        DO sta = 1,nsta
                 
          eln = elsta(sta)                 
          et = el_type(eln)                 
                 
          IF (mod(et,2) == 1) THEN
            CALL tri_basis(p,ndof(et),1,rssta(sta,1),rssta(sta,2),phi)       
          ELSE IF (mod(et,2) == 0) THEN
            CALL quad_basis(p,ndof(et),1,rssta(sta,1),rssta(sta,2),phi)
          ENDIF          
          
          Hsta = 0d0
          Qxsta = 0d0
          Qysta = 0d0
          DO dof = 1,ndof(et)
            Hsta = Hsta + H(eln,dof)*phi(dof)
            Qxsta = Qxsta + Qx(eln,dof)*phi(dof)            
            Qysta = Qysta + Qy(eln,dof)*phi(dof)
          ENDDO                
          
          WRITE(61,*) xysta(1,sta),xysta(2,sta),Hsta,hbsta(sta)
          WRITE(621,*) xysta(1,sta),xysta(2,sta),Qxsta
          WRITE(622,*) xysta(1,sta),xysta(2,sta),Qysta
        ENDDO
        
      ENDDO
      
      CLOSE(63)
      CLOSE(641)
      CLOSE(642)
      
      CLOSE(61)
      CLOSE(621)
      CLOSE(622)      
      END PROGRAM stations