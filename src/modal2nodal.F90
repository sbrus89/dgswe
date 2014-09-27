      SUBROUTINE modal2nodal()

      USE globals, ONLY: pres,nel_type,p,np,nnds,mnnds,ndof,mndof,m2n,out_direc, &
                         nqpta,mnqpta,wpta,phia,phil
      USE basis, ONLY: tri_basis,quad_basis,tri_nodes,quad_nodes,linear

      IMPLICIT NONE
      
      INTEGER :: i,j,m,dof
      INTEGER :: et,pt
      INTEGER :: n
      INTEGER :: alloc_status
      REAL(pres) :: r(9),s(9)
      REAL(pres) :: phi(9*mndof)
      
      REAL(pres) :: ml2(3,ndof(1)),mml(3,3)
      REAL(pres) :: qint      

      
      ALLOCATE(m2n(mnnds,mndof,nel_type))      
      
      OPEN(unit=111,file=trim(out_direc) // "modal2nodal.d")
      
      DO et = 1,nel_type
        
        n = nnds(et)
        
        IF (mod(et,2) == 1) THEN       
        
          CALL tri_nodes(0,np(et),n,r,s)  ! get vertex coordinates
          CALL tri_basis(p,ndof(et),n,r,s,phi)              

        ELSE IF (mod(et,2) == 0) THEN

          CALL quad_nodes(0,np(et),n,r,s)  ! get vertex coordinates
          CALL quad_basis(p,ndof(et),n,r,s,phi)                             
          
        ENDIF
        
        DO dof = 1,ndof(et)            
          DO pt = 1,n 
            i = (dof-1)*n + pt   
            IF (phi(i) < 1d-17) THEN
              m2n(pt,dof,et) = 0d0
            ELSE   
              m2n(pt,dof,et) = phi(i)
            ENDIF
          ENDDO
        ENDDO
        
        WRITE(111,"(I5,I5)") n, ndof(et)
        PRINT*, n
        DO pt = 1,n
          WRITE(111,"(160(e25.17,1x))") (m2n(pt,dof,et), dof = 1,ndof(et))
        ENDDO        
          
      ENDDO
       
      CLOSE(111)
        

        
      ! Write out L2 projection information (for use with pdeplot, only for triangles)        
      
        
      ! Calculate linear nodal basis functions
      CALL linear(phil)

      ! Compute RHS L2 projection matrix
      DO i = 1,3
        DO j = 1,ndof(1)
          qint = 0d0
          DO pt = 1,nqpta(1)
            qint = qint + wpta(pt,1)*phil(i,pt,1)*phia(j,pt,1)
          ENDDO
          ml2(i,j) = qint
        ENDDO
      ENDDO

      ! Compute linear mass matrix
      DO i = 1,3
        DO j = 1,3
          qint = 0d0
          DO pt = 1,nqpta(1)
            qint = qint + wpta(pt,1)*phil(i,pt,1)*phil(j,pt,1)
          ENDDO
          mml(i,j) = qint
        ENDDO
      ENDDO


      OPEN(unit=10,file=trim(out_direc) // 'projection.d')

      WRITE(10,*) ndof(1)
      DO i = 1,3
        WRITE(10,"(160(e24.17,1x))") (ml2(i,j),j=1,ndof(1))
      ENDDO

      DO i = 1,3
        WRITE(10,"(3(e24.17,1x))") (mml(i,j),j=1,3)
      ENDDO        
            
      CLOSE(10)

      RETURN
      END SUBROUTINE modal2nodal