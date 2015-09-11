      SUBROUTINE boundary_edge_elev()

      USE globals, ONLY: pres,el_type,ndof,nqpte, &
                         nobed,obedn,nobfr, &
                         tstage,ramp, &
                         obfreq,obper,obeq, &
                         obnfact,obamp_qpt,obph_qpt, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Zqpt,Qxqpt,Qyqpt, &
                         Zflux,Qxflux,Qyflux

      IMPLICIT NONE

      INTEGER :: ed,pt,bfr
      INTEGER :: ged
      INTEGER :: led_in,el_in,gp_in
      REAL(pres) :: Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex
      REAL(pres) :: Zhat,Qxhat,Qyhat
      REAL(pres) :: tx,ty,nx,ny,nx2,ny2,nxny,sp,hb
      REAL(pres) :: arg,Qn,Qt
      
!$OMP do

      ! Open boundary edges (elevation specified)
      DO ed = 1,nobed
            
        ged = obedn(ed)
        led_in = ged2led(1,ged)
        el_in = gel2ael(ged2el(1,ged))
              
        DO pt = 1,nqpte(1)

          gp_in = (led_in-1)*nqpte(1) + pt

          nx = nx_pt(ged,pt)
          ny = ny_pt(ged,pt)
          sp = Spe(ged,pt)
              
          hb = hbqpted(ged,pt)                 

          Z_in = Zqpt(el_in,gp_in)

          Z_ex = 0d0
          DO bfr = 1,nobfr
            arg = obfreq(bfr)*(tstage-real(INT(tstage/obper(bfr)),pres)*obper(bfr)) + obeq(bfr)
            Z_ex = Z_ex + obamp_qpt(bfr,pt,ed)*obnfact(bfr)*ramp*COS(arg-obph_qpt(bfr,pt,ed))
          ENDDO

          Qx_in = Qxqpt(el_in,gp_in)
          Qy_in = Qyqpt(el_in,gp_in)

          Qx_ex = Qx_in
          Qy_ex = Qy_in
 
          CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)

          Zflux(el_in,gp_in) = detJe(ged,pt)*Zhat
          Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat
          Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat
                          
          ENDDO
        ENDDO
!$OMP end do  

      END SUBROUTINE boundary_edge_elev