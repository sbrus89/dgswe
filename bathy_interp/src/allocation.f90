      MODULE allocation

      USE globals, ONLY: grid,nverts
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes(mesh)
      
      IMPLICIT NONE   
      
      INTEGER :: hbp,ctp
      TYPE(grid) :: mesh
      
      ctp = mesh%ctp      
      hbp = mesh%hbp

      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4          
      
      mesh%np(1) = 1
      mesh%np(2) = 1
      mesh%np(3) = ctp
      mesh%np(4) = ctp   
      mesh%np(5) = hbp
      mesh%np(6) = hbp       
      
      mesh%nnds(1) = 3
      mesh%nnds(2) = 4
      mesh%nnds(3) = (ctp+1)*(ctp+2)/2
      mesh%nnds(4) = (ctp+1)*(ctp+1) 
      mesh%nnds(5) = (hbp+1)*(hbp+2)/2
      mesh%nnds(6) = (hbp+1)*(hbp+1)      
      mesh%mnnds = maxval(mesh%nnds)         

      RETURN
      END SUBROUTINE sizes      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      



      END MODULE allocation

