      MODULE evaluate
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      
      SUBROUTINE vandermonde()
      
      USE globals, ONLY: pres,ctp,np,nnds,mnnds,nel_type,V,ipiv
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n
      REAL(pres) :: r(mnnds),s(mnnds)
      REAL(pres) :: phi(mnnds*mnnds)
      INTEGER :: info      
      
      ALLOCATE(V(mnnds,mnnds,nel_type))
      ALLOCATE(ipiv(mnnds,nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(et),n,r,s)
          CALL tri_basis(np(et),n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(et),n,r,s)
          CALL quad_basis(np(et),n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            V(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE normals()
      
      USE globals, ONLY: pres,ne,el_type,np,nnds,mnnds,V,l,elxy,elhb,nhb,xc,yc,hc,ipiv
      USE basis, ONLY: tri_basis,quad_basis      

      IMPLICIT NONE
      INTEGER :: it,el,et,p,n,i,np1
      INTEGER :: info
      REAL(pres) :: x,y
      REAL(pres) :: r(1),s(1),hb
      REAL(pres) :: dxdr,dxds,dydr,dyds,jac
      REAL(pres) :: drdx,drdy,dsdx,dsdy
      REAL(pres) :: dhdx,dhdy
      REAL(pres) :: phi(mnnds),dpdr(mnnds),dpds(mnnds)
      REAL(pres) :: nx,ny,nz,nrm
      REAL(pres) :: x1,x2,y1,y2,z1,z2
      
      info = 0
      
      ALLOCATE(nhb(3,ne))
      ALLOCATE(xc(ne),yc(ne),hc(ne))
      
      DO el = 1,ne
      
      et = el_type(el)
      p = np(et)  
      n = nnds(et)
      np1 = n+1             
        
      dxdr = 0d0
      dxds = 0d0
      dydr = 0d0
      dyds = 0d0
      xc(el) = 0d0
      yc(el) = 0d0        
        
      DO i = 1,n         
        dxdr = dxdr + l(i,2*np1,et)*elxy(i,el,1)
        dxds = dxds + l(i,3*np1,et)*elxy(i,el,1)
        dydr = dydr + l(i,2*np1,et)*elxy(i,el,2)
        dyds = dyds + l(i,3*np1,et)*elxy(i,el,2)
          
        xc(el) = xc(el) + l(i,np1,et)*elxy(i,el,1)
        yc(el) = yc(el) + l(i,np1,et)*elxy(i,el,2)
      ENDDO
        
      jac = dxdr*dyds - dydr*dxds

      drdx =  dyds/jac
      drdy = -dxds/jac
      dsdx = -dydr/jac
      dsdy =  dxdr/jac
      
      dhdx = 0d0
      dhdy = 0d0
      
      hc(el) = 0d0
      DO i = 1,n
        dhdx = dhdx + (l(i,2*np1,et)*drdx + l(i,3*np1,et)*dsdx)*elhb(i,el)
        dhdy = dhdy + (l(i,2*np1,et)*drdy + l(i,3*np1,et)*dsdy)*elhb(i,el)
        
        hc(el) = hc(el) + l(i,np1,et)*elhb(i,el)
      ENDDO
      
      nrm = sqrt(dhdx**2 + dhdy**2 + 1d0)
      
      nhb(1,el) = dhdx/nrm
      nhb(2,el) = dhdy/nrm
      nhb(3,el) = -1d0/nrm
      
      x1 = elxy(2,el,1) - elxy(1,el,1)
      x2 = elxy(3,el,1) - elxy(1,el,1)
      
      y1 = elxy(2,el,2) - elxy(1,el,2)
      y2 = elxy(3,el,2) - elxy(1,el,2)
      
      z1 = elhb(2,el) - elhb(1,el)
      z2 = elhb(3,el) - elhb(1,el)
      
      nx =   y1*z2 - y2*z1
      ny = -(x1*z2 - x2*z1)
      nz =   x1*y2 - x2*y1
      
      nrm = sqrt(nx**2 + ny**2 + nz**2)
      
      nx = nx/nrm
      ny = ny/nrm
      nz = nz/nrm
      
      PRINT*, (nhb(i,el), i=1,3)
      PRINT*, nx,ny,nz
      PRINT*, " " 
        
      ENDDO  
      
      RETURN 
      END SUBROUTINE normals
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE coordinates()
      
      USE globals, ONLY: pres,ne,mnnds,nnds,el_type,l,elxy,elhb,xpt,ypt,hpt
      
      IMPLICIT NONE
      
      INTEGER :: el,i,pt
      INTEGER :: et,n,p
      
      ALLOCATE(xpt(mnnds,ne),ypt(mnnds,ne),hpt(mnnds,ne))
      
      DO el = 1,ne
      
        et = el_type(el)
        IF (mod(et,2) == 1) THEN
          et = 3
          n = nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          et = 4
          n = nnds(4)
        ENDIF            
        
        DO pt = 1,n
          xpt(pt,el) = 0d0
          ypt(pt,el) = 0d0
          hpt(pt,el) = 0d0
          DO i = 1,n                   
            xpt(pt,el) = xpt(pt,el) + l(i,pt,et)*elxy(i,el,1)
            ypt(pt,el) = ypt(pt,el) + l(i,pt,et)*elxy(i,el,2)
            hpt(pt,el) = hpt(pt,el) + l(i,pt,et)*elhb(i,el)
          ENDDO        
        ENDDO
      
      ENDDO
      
      RETURN
      END SUBROUTINE coordinates
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE transformation()
      
      USE globals, ONLY: pres,mnnds,nel_type,nnds,np,l,V,ipiv
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: et,pt,i,j
      INTEGER :: n,p,np1
      INTEGER :: info
      REAL(pres), DIMENSION(mnnds+1,nel_type) :: r,s
      REAL(pres), DIMENSION(mnnds*(mnnds+1)) :: phi,dpdr,dpds
      
      ALLOCATE(l(mnnds,3*(mnnds+1),nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        p = np(et)
      
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,p,n,r(1,et),s(1,et))
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,p,n,r(1,et),s(1,et))
        ENDIF     
        
        np1 = n+1
        
        IF (mod(et,2) == 1) THEN
          r(np1,et) = -1d0/3d0
          s(np1,et) = -1d0/3d0
        ELSE IF (mod(et,2) == 0) THEN
          r(np1,et) = 1d0
          s(np1,et) = 1d0
        ENDIF        
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,np1,r(1,et),s(1,et),phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,np1,r(1,et),s(1,et),phi,dpdr,dpds)
        ENDIF     
        
        DO pt = 1,np1
          DO i = 1,n
            j = (i-1)*np1+pt
          
            l(i,pt,et) = phi(j)
            l(i,np1+pt,et) = dpdr(j)
            l(i,2*np1+pt,et) = dpds(j)          
          ENDDO         
        ENDDO

        CALL DGETRS("N",n,3*np1,V(1,1,et),mnnds,ipiv(1,et),l(1,1,et),mnnds,info)
      ENDDO
      
      
      RETURN
      END SUBROUTINE transformation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      
      END MODULE evaluate