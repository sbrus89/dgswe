      SUBROUTINE curvilinear(mesh)

      USE globals, ONLY: rp,nel_type,np,ctp,V,ipiv,nnds,nverts,mnnds,grid
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis

      IMPLICIT NONE
      
      TYPE(grid) :: mesh      
      
      INTEGER :: ed,led,nd,pt
      INTEGER :: ged,el,ind
      INTEGER :: et,etl
      INTEGER :: nvert,nnd,npts,n
      INTEGER :: p
      INTEGER :: info
      REAL(rp) :: xpt,ypt,ytest,hb
      REAL(rp) :: r(mnnds),s(mnnds)        
      REAL(rp) :: x(mnnds),y(mnnds)      
      REAL(rp) :: l(mnnds,mnnds,nel_type)
      

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluate linear shape functions at curvilinear nodal sets (to compute curvilnear nodes)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      DO et = 1,nel_type      
        
        p = np(et)
        
        IF (mod(et,2) == 1) THEN
          etl = 1
          CALL tri_nodes(1,p,npts,r,s)
          CALL tri_basis(np(etl),n,npts,r,s,l(:,:,et))
        ELSE IF (mod(et,2) == 0) THEN
          etl = 2
          CALL quad_nodes(1,p,npts,r,s)
          CALL quad_basis(np(etl),n,npts,r,s,l(:,:,et))
        ENDIF       
        
        CALL DGETRS("N",n,npts,V(1,1,etl),mnnds,ipiv(1,etl),l(1,1,et),mnnds,info)            
      
      ENDDO      
      
      
      
      
      
      DO ed = 1,mesh%nnfbed
        ged = mesh%nfbedn(ed)
        
        el = mesh%ged2el(1,ged)
        led = mesh%ged2led(1,ged)
        
        IF (mesh%nelnds(el) == 3) THEN
          mesh%el_type(el) = 3 
          etl = 1
        ELSE IF (mesh%nelnds(el) == 4) THEN
          mesh%el_type(el) = 4        
          etl = 2
        ENDIF
        
        et = mesh%el_type(el)        
        nnd = nnds(et)
        nvert = nverts(et)
        
        mesh%nelnds(el) = nnd         
        
        DO pt = 1,nnd     
          
          DO nd = 1,nvert
            x(nd) = mesh%xy(1,mesh%ect(nd,el))
            y(nd) = mesh%xy(2,mesh%ect(nd,el))
          ENDDO
          
          xpt = 0d0
          ypt = 0d0
          
          DO nd = 1,nvert
            xpt = xpt + l(nd,pt,et)*x(nd)
            ypt = ypt + l(nd,pt,et)*y(nd)
          ENDDO
          
          mesh%elxy(pt,el,1) = xpt
          mesh%elxy(pt,el,2) = ypt
        ENDDO      
        
        
        
        DO nd = 1,ctp-1
          pt = mod(led,nvert)*ctp + 1 + nd
          
          ytest = mesh%elxy(pt,el,2)
          xpt = mesh%elxy(pt,el,1) 
          
          IF (ytest < 250d0) THEN
            ypt = 0d0 + 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
          ELSE IF (ytest > 250d0) THEN
            ypt = 500d0 - 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
          ENDIF
          
          mesh%elxy(pt,el,2) = ypt
        ENDDO       
        
        
        
      
        
        
      ENDDO
      

      
      
      RETURN
      END SUBROUTINE curvilinear
