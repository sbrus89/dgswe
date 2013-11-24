      SUBROUTINE numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)
      
!       USE globals, ONLY:pres,nx,ny,g,eig_in,eig_ex,alpha, &
!                         H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!                         xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!                         Hhat,Qxhat,Qyhat

      USE globals, ONLY:pres,g


      IMPLICIT NONE
      
      REAL(pres) :: eig_in,eig_ex,alpha
      REAL(pres) :: Qx_in,Qy_in,H_in
      REAL(pres) :: Qx_ex,Qy_ex,H_ex
      REAL(pres) :: nx,ny
      REAL(pres) :: Qxhat,Qyhat,Hhat
      REAL(pres) :: xmom_in,ymom_in,xymom_in
      REAL(pres) :: xmom_ex,ymom_ex,xymom_ex
      REAL(pres) :: rHin,rHex
      
      rHin = 1d0/H_in
      rHex = 1d0/H_ex
      
      eig_in = abs(Qx_in*nx + Qy_in*ny)*rHin + sqrt(g*H_in) 
      eig_ex = abs(Qx_ex*nx + Qy_ex*ny)*rHex + sqrt(g*H_ex)

      alpha = max(eig_in,eig_ex)

!       eig(1) = abs((Qx_in*nx + Qy_in*nx)/H_in + sqrt(g*H_in))
!       eig(2) = abs((Qx_in*nx + Qy_in*nx)/H_in)
!       eig(3) = abs((Qx_in*nx + Qy_in*nx)/H_in - sqrt(g*H_in))
! 
!       eig(4) = abs((Qx_ex*nx + Qy_ex*nx)/H_ex + sqrt(g*H_ex))
!       eig(5) = abs((Qx_ex*nx + Qy_ex*nx)/H_ex)
!       eig(6) = abs((Qx_ex*nx + Qy_ex*nx)/H_ex - sqrt(g*H_ex))
! 
!       alpha = maxval(eig,1)

      xmom_in = Qx_in*Qx_in*rHin + .5d0*g*H_in*H_in
      xmom_ex = Qx_ex*Qx_ex*rHex + .5d0*g*H_ex*H_ex
      
      ymom_in = Qy_in*Qy_in*rHin + .5d0*g*H_in*H_in
      ymom_ex = Qy_ex*Qy_ex*rHex + .5d0*g*H_ex*H_ex
       
      xymom_in = Qx_in*Qy_in*rHin
      xymom_ex = Qx_ex*Qy_ex*rHex

      Hhat = .5d0*( nx*(Qx_in+Qx_ex) + ny*(Qy_in+Qy_ex) - alpha*(H_ex-H_in) )
      Qxhat = .5d0*( nx*(xmom_in+xmom_ex) + ny*(xymom_in+xymom_ex) - alpha*(Qx_ex-Qx_in) )
      Qyhat = .5d0*( nx*(xymom_in+xymom_ex) + ny*(ymom_in+ymom_ex) - alpha*(Qy_ex-Qy_in) )

      RETURN
      END SUBROUTINE numerical_flux