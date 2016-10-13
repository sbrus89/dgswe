      MODULE evaluate
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE function_eval()
      
      USE globals, ONLY: rp,nel_type,nqpta,qpta,mnqpta,fine
      USE basis, ONLY: tri_basis,quad_basis,dgswem_basis,linear_basis
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      
      IMPLICIT NONE
      
      INTEGER :: p,n,et,m,i,pt,mnnds,mndof,npt
      INTEGER :: info
      REAL(rp) :: r(mnqpta),s(mnqpta)
      
      mnnds = fine%mnnds
      mndof = fine%mndof         
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluate shape functions at quadrature points (to compute r,s -> x,y transformtion in error integration)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(fine%l(mnnds,mnqpta,nel_type))   
      ALLOCATE(fine%dldr(mnnds,mnqpta,nel_type))
      ALLOCATE(fine%dlds(mnnds,mnqpta,nel_type))
      
      DO et = 1,nel_type
        p = fine%np(et)     ! transformation order
        npt = nqpta(et)
        
        DO pt = 1,npt
          r(pt) = qpta(pt,1,et) 
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        CALL shape_functions_area_eval(et,p,n,npt,r,s,fine%l(:,:,et),fine%dldr(:,:,et),fine%dlds(:,:,et))             
      
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluate basis functions at quadrature points (to evaluate solution in error integration)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(fine%phi(mndof,mnqpta,nel_type))  
      
      DO et = 1,nel_type
        p = fine%p        ! solution order
        npt = nqpta(et) 
        
        DO pt = 1,npt
          r(pt) = qpta(pt,1,et) 
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
#ifdef dgswem
          CALL dgswem_basis(p,n,npt,r,s,fine%phi(:,:,et))
#elif adcirc           
          CALL linear_basis(npt,r,s,fine%phi(:,:,et))   
#else
          CALL tri_basis(p,n,npt,r,s,fine%phi(:,:,et))          
#endif          
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,fine%phi(:,:,et))
        ENDIF        
      
      ENDDO      
      
      RETURN
      END SUBROUTINE function_eval
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE detJ_eval()
      
      USE globals, ONLY: rp,fine,mnqpta,nqpta
      USE transformation, ONLY: element_transformation
      
      IMPLICIT NONE
      
      INTEGER :: el,pt,nd
      INTEGER :: et
      REAL(rp) :: xpt,ypt
      REAL(rp) :: drdx,dsdx,drdy,dsdy      
      
      ALLOCATE(fine%detJ(mnqpta,fine%ne))
      
      ! Calculate the determinant of the Jacobian at quadrature points (to compute integal in error integration)
      DO el = 1,fine%ne
        et = fine%el_type(el)
        
        DO pt = 1,nqpta(et)
    
          
          CALL element_transformation(fine%nnds(et),fine%elxy(:,el,1),fine%elxy(:,el,2),fine%l(:,pt,et),xpt,ypt, &
                                      fine%dldr(:,pt,et),fine%dlds(:,pt,et),drdx,drdy,dsdx,dsdy,fine%detJ(pt,el))
          
        ENDDO
      
      ENDDO
      
      RETURN
      END SUBROUTINE detJ_eval      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      END MODULE evaluate