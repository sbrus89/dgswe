      MODULE allocation

      USE globals, ONLY: solution,nverts
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes(sol)
      
      IMPLICIT NONE   
      
      INTEGER :: p,ctp
      TYPE(solution) :: sol
      
      p = sol%p
      ctp = sol%ctp

      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      
      
      sol%ndof(1) = (p+1)*(p+2)/2
      sol%ndof(2) = (p+1)*(p+1)
      sol%ndof(3) = sol%ndof(1)
      sol%ndof(4) = sol%ndof(2)
      sol%mndof = maxval(sol%ndof)
      
      sol%np(1) = 1
      sol%np(2) = 1
      sol%np(3) = ctp
      sol%np(4) = ctp        
      
      sol%nnds(1) = 3
      sol%nnds(2) = 4
      sol%nnds(3) = (ctp+1)*(ctp+2)/2
      sol%nnds(4) = (ctp+1)*(ctp+1) 
      sol%mnnds = maxval(sol%nnds)         

      RETURN
      END SUBROUTINE sizes      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE nest_alloc()
      
      USE globals, ONLY: fine,elf2elc,elf2elb 
      
      IMPLICIT NONE
      
      ALLOCATE(elf2elc(fine%ne),elf2elb(fine%ne))    
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "          Finding Element Nesting            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "        
      
      RETURN
      END SUBROUTINE nest_alloc
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      SUBROUTINE error_alloc()
      
      USE globals, ONLY: mnqpta,coarse,fine, &
                         Hc,Qxc,Qyc,Hf,Qxf,Qyf, &
                         xf,yf,r,s,hb,phi
      
      IMPLICIT NONE
      
      INTEGER :: ne,mndof
      
      ALLOCATE(xf(mnqpta),yf(mnqpta))
      ALLOCATE(r(mnqpta),s(mnqpta),hb(mnqpta))
      
      mndof = coarse%mndof
      ne = coarse%ne
      
!       ALLOCATE(coarse%H(ne,mndof,1),coarse%Qx(ne,mndof,1),coarse%Qy(ne,mndof,1))
      ALLOCATE(coarse%phi(mndof,mnqpta,1))
      ALLOCATE(phi(mndof*mnqpta))
 
      mndof = fine%mndof
      ne = fine%ne
      
!       ALLOCATE(fine%H(ne,mndof,1),fine%Qx(ne,mndof,1),fine%Qy(ne,mndof,1)) 
      
      ALLOCATE(Hc(mnqpta),Qxc(mnqpta),Qyc(mnqpta))
      ALLOCATE(Hf(mnqpta),Qxf(mnqpta),Qyf(mnqpta))    
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Reading Solutions                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      RETURN
      END SUBROUTINE  error_alloc

      END MODULE allocation

