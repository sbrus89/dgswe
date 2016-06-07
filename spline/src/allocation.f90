      MODULE allocation

      USE globals, ONLY: grid

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes(mesh)
      
      USE globals, ONLY: nverts
      USE basis, ONLY: lglpts           
      
      IMPLICIT NONE   
      
      TYPE(grid) :: mesh
      
      mesh%np(1) = 1
      mesh%np(2) = 1
      mesh%np(3) = mesh%ctp
      mesh%np(4) = mesh%ctp        
      
      mesh%nnds(1) = 3
      mesh%nnds(2) = 4
      mesh%nnds(3) = (mesh%ctp+1)*(mesh%ctp+2)/2
      mesh%nnds(4) = (mesh%ctp+1)*(mesh%ctp+1) 
      mesh%mnnds = maxval(mesh%nnds)    
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      

      ALLOCATE(mesh%rpts(mesh%ctp+1))
      
!       DO j = 0,mesh%ctp
!         mesh%rpts(j) = -1d0 + real(j,rp)*2d0/real(mesh%ctp,rp)   
!       ENDDO
      
      CALL lglpts(mesh%ctp,mesh%rpts)
            
      
      RETURN
      END SUBROUTINE sizes      
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
  
      END MODULE allocation