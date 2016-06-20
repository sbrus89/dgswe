      SUBROUTINE sizes(ndof,mndof,nverts,np,mnp,nnds,mnnds)
                             
      USE read_dginp, ONLY: p,ctp,hbp
      
      IMPLICIT NONE
      
      INTEGER :: ndof(4)
      INTEGER :: mndof
      INTEGER :: nverts(4)
      INTEGER :: np(4)
      INTEGER :: mnp
      INTEGER :: nnds(4)
      INTEGER :: mnnds
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)      
      mndof = maxval(ndof)
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp  
      mnp = maxval(np)+1

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      mnnds = maxval(nnds)      
      
            
      RETURN
      END SUBROUTINE sizes