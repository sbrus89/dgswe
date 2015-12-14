      MODULE evaluate
      
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE vandermonde(mesh)
      
      USE globals, ONLY: rp,nel_type,grid
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      TYPE(grid) :: mesh
      
      INTEGER :: et,n
      REAL(rp) :: r(mesh%mnnds),s(mesh%mnnds)
      INTEGER :: info    
      

      
      ALLOCATE(mesh%V(mesh%mnnds,mesh%mnnds,nel_type+2))
      ALLOCATE(mesh%ipiv(mesh%mnnds,nel_type+2))
      
      ! Evaluate basis functions at reference element nodes
      DO et = 1,nel_type+2
        n = mesh%nnds(et)

        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,mesh%np(et),n,r,s)
          CALL tri_basis(mesh%np(et),n,n,r,s,mesh%V(:,:,et))       
        ELSE 
          CALL quad_nodes(1,mesh%np(et),n,r,s)
          CALL quad_basis(mesh%np(et),n,n,r,s,mesh%V(:,:,et))
        ENDIF
        
        
        ! Do LU decomposition 
        CALL DGETRF(n,n,mesh%V(1,1,et),mesh%mnnds,mesh%ipiv(1,et),info)
        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (mesh%V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "       
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE re_vert(mesh)
      
      USE globals, ONLY: rp,grid,nel_type
      USE basis, ONLY: tri_nodes, quad_nodes
      
      IMPLICIT NONE
      
      INTEGER :: et,pt
      INTEGER :: n,p
      REAL(rp) :: r(4),s(4)
      TYPE(grid) :: mesh
      
      DO et = 1,nel_type
        IF (mod(et,2) == 1) THEN
          n = mesh%nnds(1)
          p = mesh%np(1)
          CALL tri_nodes(1,p,n,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          n = mesh%nnds(2)
          p = mesh%np(2)
          CALL quad_nodes(1,p,n,r,s)
        ENDIF
        
        DO pt = 1,n
          mesh%rsre(1,pt,et) = r(pt)
          mesh%rsre(2,pt,et) = s(pt)
        ENDDO
        
!         PRINT("(4(F15.5))"), (mesh%rsre(1,pt,et), pt = 1,n)
!         PRINT("(4(F15.5))"), (mesh%rsre(2,pt,et), pt = 1,n)        
!         PRINT*, " "
        
      ENDDO      
      
      RETURN
      END SUBROUTINE re_vert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE function_eval()
      
      USE globals, ONLY: rp,nel_type,mnept,nept,ept,eval,base
      USE basis, ONLY: tri_basis,quad_basis,adcirc_basis
      
      IMPLICIT NONE
           
      INTEGER :: p,n,et,m,i,pt,mnnds,mndof,npt
      INTEGER :: info
      
      REAL(rp) :: r(mnept),s(mnept)
      
      mnnds = eval%mnnds
          
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluate shape functions at evaluation points (to compute r,s -> x,y transformtion)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(eval%l(mnnds,mnept,nel_type))   
      ALLOCATE(eval%dldr(mnnds,mnept,nel_type))
      ALLOCATE(eval%dlds(mnnds,mnept,nel_type))
      
      DO et = 1,nel_type
        p = eval%np(et)     ! transformation order
        n = eval%nnds(et)   ! transformation nodes
        npt = nept(et)
        
        DO pt = 1,npt
          r(pt) = ept(pt,1,et) 
          s(pt) = ept(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,eval%l(:,:,et),eval%dldr(:,:,et),eval%dlds(:,:,et))
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,eval%l(:,:,et),eval%dldr(:,:,et),eval%dlds(:,:,et))
        ENDIF       
        
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%l(1,1,et),mnnds,info)  
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%dldr(1,1,et),mnnds,info)      
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%dlds(1,1,et),mnnds,info)              
      
      ENDDO
      
      


      
      RETURN
      END SUBROUTINE function_eval
      
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
 
      SUBROUTINE eval_hb(ele,pte,xe,x,hb)

      USE globals, ONLY: rp,base,eval
      USE basis, ONLY: tri_basis,quad_basis
      USE find_element, ONLY: in_element

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ele
      INTEGER, INTENT(IN) :: pte
      REAL(rp), INTENT(IN), OPTIONAL :: xe(2)      
      INTEGER :: elb,ete,etb,et
      INTEGER :: npts,nd,n
      INTEGER :: info
      INTEGER :: exceed
      REAL(rp), INTENT(OUT) :: x(2)
      REAL(rp), INTENT(OUT) :: hb      
      REAL(rp) :: r(2)
      REAL(rp) :: error

      ete = eval%el_type(ele)
      
      IF (PRESENT(xe)) THEN
        ! Take x,y coordinates as evaluation point (ele and pte are unused)
        x(1) = xe(1)      
        x(2) = xe(2)
      ELSE
        ! evaluate x,y coordinates of eval element evaluation point
        x(1) = 0d0
        x(2) = 0d0
      
        DO nd = 1,eval%nnds(ete)
          x(1) = x(1) + eval%l(nd,pte,ete)*eval%elxy(nd,ele,1)  
          x(2) = x(2) + eval%l(nd,pte,ete)*eval%elxy(nd,ele,2)
        ENDDO   
      ENDIF
          
!       PRINT*,ele,pte
!       PRINT*,x(1),x(2)
      CALL in_element(x(1:2),elb,r(1:2),error,exceed)     
          
      etb = base%el_type(elb)  ! element type for base element          
      
      IF (exceed == 1) THEN
        PRINT*, "MAX ITERATIONS EXCEEDED IN FINDING R,S COORDINATES"
        PRINT*, "ERROR: ",error
        PRINT*, "ELEMENT TYPE: ", etb
      ELSE 
!         PRINT*, "error: ",error
!         PRINT*, "element type: ", etb        
      ENDIF
      
        
      IF (mod(etb,2) == 1) THEN
        et = 5
        CALL tri_basis(base%hbp,n,1,r(1),r(2),base%l(:,:,et))
      ELSE IF (mod(etb,2) == 0) THEN
        et = 6
        CALL quad_basis(base%hbp,n,1,r(1),r(2),base%l(:,:,et))
      ENDIF       
        
      CALL DGETRS("N",n,1,base%V(1,1,et),base%mnnds,base%ipiv(1,et),base%l(1,1,et),base%mnnds,info)                          
               
               
      ! Evaluate bathymetry at r,s coordinates
      hb = 0d0
      DO nd = 1,n
        hb = hb + base%l(nd,1,et)*base%elhb(nd,elb)
      ENDDO            

  

      RETURN
      END SUBROUTINE eval_hb


      
      
      END MODULE evaluate