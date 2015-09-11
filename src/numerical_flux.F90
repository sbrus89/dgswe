      SUBROUTINE numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
      
!       USE globals, ONLY:rp,nx,ny,g,eig_in,eig_ex,alpha, &
!                         H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!                         xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!                         Hhat,Qxhat,Qyhat

      USE globals, ONLY:rp,g


      IMPLICIT NONE
      
      REAL(rp) :: eig_in,eig_ex,alpha
      REAL(rp) :: Qx_in,Qy_in,Z_in
      REAL(rp) :: Qx_ex,Qy_ex,Z_ex
      REAL(rp) :: nx,ny,sp,hb
      REAL(rp) :: Qxhat,Qyhat,Zhat
      REAL(rp) :: H_in,H_ex,hb_in,hb_ex
      REAL(rp) :: xmom_in,ymom_in,xymom_in
      REAL(rp) :: xmom_ex,ymom_ex,xymom_ex
      REAL(rp) :: rHin,rHex
      REAL(rp) :: nxsp,cfac
      
      hb_in = hb
      hb_ex = hb
      
      H_in = Z_in + hb_in
      H_ex = Z_ex + hb_ex
            
      rHin = 1d0/H_in
      rHex = 1d0/H_ex
      
      nxsp = nx*sp
      cfac = ny**2 + (nxsp)**2
      
      eig_in = abs(Qx_in*nxsp + Qy_in*ny)*rHin + sqrt(g*H_in*cfac) 
      eig_ex = abs(Qx_ex*nxsp + Qy_ex*ny)*rHex + sqrt(g*H_ex*cfac)

      alpha = max(eig_in,eig_ex)

!       eig(1) = abs((Qx_in*nxsp + Qy_in*ny)/H_in + sqrt(g*H_in*cfac))
!       eig(2) = abs((Qx_in*nxsp + Qy_in*ny)/H_in)
!       eig(3) = abs((Qx_in*nxsp + Qy_in*ny)/H_in - sqrt(g*H_in*cfac))
! 
!       eig(4) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex + sqrt(g*H_ex*cfac))
!       eig(5) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex)
!       eig(6) = abs((Qx_ex*nxsp + Qy_ex*ny)/H_ex - sqrt(g*H_ex*cfac))
! 
!       alpha = maxval(eig,1)

      xmom_in = Qx_in*Qx_in*rHin + .5d0*g*(H_in*H_in - hb_in*hb_in)
      xmom_ex = Qx_ex*Qx_ex*rHex + .5d0*g*(H_ex*H_ex - hb_ex*hb_ex)
      
      ymom_in = Qy_in*Qy_in*rHin + .5d0*g*(H_in*H_in - hb_in*hb_in)
      ymom_ex = Qy_ex*Qy_ex*rHex + .5d0*g*(H_ex*H_ex - hb_ex*hb_ex)
       
      xymom_in = Qx_in*Qy_in*rHin
      xymom_ex = Qx_ex*Qy_ex*rHex

      Zhat  = .5d0*( nxsp*(Qx_in+Qx_ex) + ny*(Qy_in+Qy_ex) - alpha*(Z_ex-Z_in) )
      Qxhat = .5d0*( nxsp*(xmom_in+xmom_ex) + ny*(xymom_in+xymom_ex) - alpha*(Qx_ex-Qx_in) )
      Qyhat = .5d0*( nxsp*(xymom_in+xymom_ex) + ny*(ymom_in+ymom_ex) - alpha*(Qy_ex-Qy_in) )

      RETURN
      END SUBROUTINE numerical_flux