      SUBROUTINE shape_functions_at_qpts(mesh)
      
      USE globals, ONLY: rp,grid,nverts,nel_type
                         
      USE area_qpts_mod, ONLY: area_qpts
      USE edge_qpts_mod, ONLY: edge_qpts
      USE shape_functions_mod, ONLY: shape_functions_area_eval

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      
      INTEGER :: i,et
      INTEGER :: nv,p,nnd
      INTEGER :: npt
      INTEGER :: myrank
      INTEGER :: mnqpt,mnnds
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s

      myrank = 0
     
      p = mesh%ctp
      
      CALL area_qpts(myrank,p,mesh%ctp,nel_type,mesh%nqpta,mesh%mnqpta,mesh%wpta,mesh%qpta)    
      CALL edge_qpts(myrank,p,mesh%ctp,nel_type,mesh%nqpte,mesh%mnqpte,mesh%wpte,mesh%qpte) 
      
      mnqpt = mesh%mnqpta+4*mesh%mnqpte
      
      ALLOCATE(r(mnqpt),s(mnqpt))
      ALLOCATE(mesh%psi(mesh%mnnds,mnqpt,nel_type),mesh%dpdr(mesh%mnnds,mnqpt,nel_type),mesh%dpds(mesh%mnnds,mnqpt,nel_type))
      
      DO et = 1,nel_type
      
        nv = nverts(et)
        p = mesh%np(et)       
      
        npt = 0
        DO i = 1,mesh%nqpta(et)
          npt = npt + 1      
      
          r(npt) = mesh%qpta(i,1,et)
          s(npt) = mesh%qpta(i,2,et)       
        ENDDO
      
        DO i = 1,nv*mesh%nqpte(et)        
          npt = npt + 1      
      
          r(npt) = mesh%qpte(i,1,et)
          s(npt) = mesh%qpte(i,2,et)
        ENDDO           
      
        CALL shape_functions_area_eval(nv,p,nnd,npt,r,s,mesh%psi(:,:,et),mesh%dpdr(:,:,et),mesh%dpds(:,:,et))   
      
      ENDDO

      RETURN
      END SUBROUTINE shape_functions_at_qpts