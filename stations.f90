      PROGRAM stations

      USE globals
      USE kdtree2_module

      IMPLICIT NONE

      INTEGER :: el,nd,sta,i
      INTEGER :: n1,n2,nvert,mnepn,eln,clnd,et
      INTEGER :: srchdp = 1
      REAL(pres) :: x(4),y(4),area,sarea,tol,found 
      TYPE(kdtree2), POINTER :: tree
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults
      
      tol = 1d-10
      
      CALL read_input()
      CALL read_grid()
      CALL read_stations()
      
      
      
      
      ! Build kd-tree
      IF (ne < srchdp) THEN
        srchdp = ne
      ENDIF      
      ALLOCATE(kdresults(srchdp))            
      tree => kdtree2_create(xy, rearrange=.true., sort=.true.)
      
      
      
      
      
      ! Find elements associated with each node            
      ALLOCATE(nepn(nn))
      
      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
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
        DO nd = 1,nvert
          n1 = vct(nd,el)
          nepn(n1) = nepn(n1) + 1
          epn(nepn(n1),n1) = el
        ENDDO
      ENDDO      
      
      
      
      
      ! Find element each station is located in
      DO sta = 1,nsta
      
        ! Find node closest to station
        CALL kdtree2_n_nearest(tp=tree,qv=xysta(:,sta),nn=srchdp,results=kdresults)
        clnd = kdresults(1)%idx
        
        PRINT*, sta,ndsta(sta),kdresults(1)%idx,kdresults(1)%dis      
        
        ! Test elements associated with nearest node to see which element station is located in
        found = 0  
  elem: DO el = 1,nepn(clnd)
          eln = epn(el,clnd)
          
          et = el_type(eln)
          nvert = nverts(et)
          
          DO i = 1,nvert
            x(i) = xy(1,vct(i,eln))
            y(i) = xy(2,vct(i,eln))
          ENDDO
          
          ! Compute element area
          IF (mod(et,2) == 1) THEN
            area = .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
          ELSE IF (mod(et,2) == 0) THEN
            area = .5d0*abs((x(3)-x(1))*(y(4)-y(2))-(x(4)-x(2))*(y(3)-y(1)))
          ENDIF          
          
          ! Compute sum of sub-triangle areas
          sarea = 0d0
          DO i = 1,nvert
            n1 = mod(i+0,nvert)+1
            n2 = mod(i+1,nvert)+1
            
            x(1) = xy(1,vct(n1,eln))
            y(1) = xy(2,vct(n1,eln))
            
            x(2) = xy(1,vct(n2,eln))
            y(2) = xy(2,vct(n2,eln))
            
            x(3) = xysta(1,sta)
            y(3) = xysta(2,sta)
            
            sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
          ENDDO
          
          ! The station is in the element if the element area and sum of sub triangle are the same
          IF (abs(area - sarea) < tol) THEN
            PRINT*, "element found", eln
            elsta(sta) = eln
            found = 1
            EXIT elem                        
          ENDIF
          
        ENDDO elem
        
        IF (found == 0) THEN
          PRINT*, "ERROR: ELEMENT NOT FOUND"
        ENDIF
      
      ENDDO
      
      END PROGRAM stations