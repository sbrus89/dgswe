      SUBROUTINE filter_normals()
      
      USE globals, ONLY: base,tree_c,kdresults
      USE kdtree2_module 
      USE evaluate
      
      IMPLICIT NONE
      
      INTEGER :: j,i,it
      INTEGER :: nneigh,el
      REAL(rp) :: xi(3),pj(3),ni(3),nj(3),top(3),bottom
      REAL(rp) :: hpt,phiw,phi0
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: nrm
      
      ALLOCATE(nrm(3,base%ne))
      
      DO j = 1,base%ne
      
        pj = base%xyhc(:,j)       
      
        CALL kdtree2_n_nearest(tp=tree_c,qv=pj,nn=1,results=kdresults)
          
        el = kdresults(1)%idx           
        hpt = (r*base%h(el))**2   
          
        CALL kdtree2_r_nearest(tp=tree_c,qv=pj,r2=hpt,nfound=nneigh,nalloc=base%ne,results=kdresults)
          
        hpt = sqrt(hpt)  
        
        top = 0d0
        bottom = 0d0
        DO i = 1,nneigh
          xi  = base%xyhc(:,kdresults(i)%idx)
          ni = base%nhb(:,kdresults(i)%idx)  
          
          phi0 = (1d0-norm(pj-xi)**2/hpt**2)**4
          
          top = top + phi0*ni
          bottom = bottom + phi0
          
        ENDDO       
        
        nrm(:,j) = top/bottom
        
  iter: DO it = 1,maxptit
          top = 0d0
          bottom = 0d0
          nj = nrm(:,j)
          DO i = 1,nneigh
          
            xi  = base%xyhc(:,kdresults(i)%idx)
            ni = base%nhb(:,kdresults(i)%idx)   
                        
            phiw = (1d0-norm(pj-xi)**2/hpt**2)**4*exp(-(norm(nj-ni)/sigma_n)**2)
            
            top = top + phiw*ni
            bottom = bottom + phiw
          ENDDO
          nrm(:,j) = top/bottom
          IF (norm(nrm(:,j)-nj) < threshold) THEN
            PRINT*, "iter = ", it
            EXIT iter
          ENDIF
        ENDDO iter
!         IF (it == maxptit+1) THEN
!           PRINT*, "max iterations reached"
!         ENDIF
      ENDDO

      base%nhb = nrm
      
      RETURN
      END SUBROUTINE filter_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      