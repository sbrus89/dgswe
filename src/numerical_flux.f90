      SUBROUTINE numerical_flux()
      
      USE globals, ONLY:nx,ny,eig_in,eig_ex,eig,alpha,g, &
                        H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
                        xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                        Hhat,Qxhat,Qyhat

      IMPLICIT NONE

      eig_in = abs(Qx_in*nx + Qy_in*ny)/H_in + sqrt(g*H_in) 
      eig_ex = abs(Qx_ex*nx + Qy_ex*ny)/H_ex + sqrt(g*H_ex)

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

      Hhat = .5d0*( nx*(Qx_in+Qx_ex) + ny*(Qy_in+Qy_ex) - alpha*(H_ex-H_in) )
      Qxhat = .5d0*( nx*(xmom_in+xmom_ex) + ny*(xymom_in+xymom_ex) - alpha*(Qx_ex-Qx_in) )
      Qyhat = .5d0*( nx*(xymom_in+xymom_ex) + ny*(ymom_in+ymom_ex) - alpha*(Qy_ex-Qy_in) )

      RETURN
      END SUBROUTINE numerical_flux