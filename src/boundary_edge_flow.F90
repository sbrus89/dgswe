      SUBROUTINE boundary_edge_flow()

      USE globals, ONLY: rp,el_type,ndof,nqpte, &
                         nfbed,fbedn,nfbfr, &
                         tstage,ramp, &
                         fbfreq,fbper,fbeq, &
                         fbnfact,fbamp_qpt,fbph_qpt, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Zqpt,Qxqpt,Qyqpt, &
                         Zflux,Qxflux,Qyflux

      IMPLICIT NONE

      INTEGER :: ed,pt,bfr
      INTEGER :: ged
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex
      REAL(rp) :: Zhat,Qxhat,Qyhat
      REAL(rp) :: tx,ty,nx,ny,nx2,ny2,nxny,sp,hb
      REAL(rp) :: arg,Qn,Qt
      

!$OMP do

      ! Flow specified boundary edges
      DO ed = 1,nfbed
            
        ged = fbedn(ed)
        led_in = ged2led(1,ged)  
        el_in = gel2ael(ged2el(1,ged))      
              
        DO pt = 1,nqpte(1)

          gp_in = (led_in-1)*nqpte(1) + pt
       
          nx = nx_pt(ged,pt)
          ny = ny_pt(ged,pt)
          sp = Spe(ged,pt)      
              
          hb = hbqpted(ged,pt)              

          tx = -ny
          ty = nx

          Z_in = Zqpt(el_in,gp_in)
          Z_ex = Z_in

          Qx_in = Qxqpt(el_in,gp_in)
          Qy_in = Qyqpt(el_in,gp_in)

          Qn = 0d0
          Qt = 0d0
          DO bfr = 1,nfbfr
            arg = fbfreq(bfr)*(tstage - real(INT(tstage/fbper(bfr)),rp)*fbper(bfr)) + fbeq(bfr)
            Qn = Qn + fbamp_qpt(bfr,pt,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(bfr,pt,ed))
          ENDDO

          Qn = -Qn

          Qx_ex = ( ty*Qn - ny*Qt)/(nx*ty-ny*tx)
          Qy_ex = (-tx*Qn + nx*Qt)/(nx*ty-ny*tx)
 
          CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)

          Zflux(el_in,gp_in) = detJe(ged,pt)*Zhat              
          Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat
          Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
        ENDDO
      ENDDO
            
!$OMP end do

      END SUBROUTINE boundary_edge_flow