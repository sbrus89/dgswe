      SUBROUTINE boundary_edge_land()

      USE globals, ONLY: rp,el_type,ndof,nqpte,nnfbed,nfbedn, &
                         ged2led,ged2el,gel2ael, &
                         nx_pt,ny_pt,Spe,detJe,hbqpted, &
                         Z,Qx,Qy, &
                         Zqpt,Qxqpt,Qyqpt, &
                         Zflux,Qxflux,Qyflux, &
                         rhsZ,rhsQx,rhsQy, &
                         phie,phie_int
                         

      IMPLICIT NONE

      INTEGER :: ed,pt,l,dof
      INTEGER :: ged,et
      INTEGER :: led_in,el_in,gp_in
      REAL(rp) :: Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex
      REAL(rp) :: Zhat,Qxhat,Qyhat
      REAL(rp) :: nx,ny,nx2,ny2,nxny,sp,hb

      IF (nqpte(3) == nqpte(1)) THEN
!$OMP do     
        ! No normal flow boundary condition 
        DO ed = 1,nnfbed
            
          ged = nfbedn(ed)
          led_in = ged2led(1,ged)        
          el_in = gel2ael(ged2el(1,ged))  
              
          DO pt = 1,nqpte(3)

            gp_in = (led_in-1)*nqpte(3) + pt

            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt)
              
            hb = hbqpted(ged,pt)         

            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny

            Z_in = Zqpt(el_in,gp_in)
            Z_ex = Z_in              

            Qx_in = Qxqpt(el_in,gp_in)
            Qy_in = Qyqpt(el_in,gp_in)

            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

            CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)

            Zflux(el_in,gp_in) = detJe(ged,pt)*Zhat
            Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat
            Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
          ENDDO
        ENDDO
!$OMP end do
      ELSE
!$OMP do     
        ! No normal flow boundary condition 
        DO ed = 1,nnfbed
            
          ged = nfbedn(ed)
          led_in = ged2led(1,ged)            
          el_in = gel2ael(ged2el(1,ged))
          et = el_type(ged2el(1,ged)) 
              
          DO pt = 1,nqpte(et)

            gp_in = (led_in-1)*nqpte(et) + pt

            nx = nx_pt(ged,pt)
            ny = ny_pt(ged,pt)
            sp = Spe(ged,pt) 
                
            hb = hbqpted(ged,pt)                

            nx2 = nx*nx
            ny2 = ny*ny
            nxny = nx*ny

            Z_in = Z(el_in,1)
            Qx_in = Qx(el_in,1)
            Qy_in = Qy(el_in,1)
            DO dof = 2,ndof(et)
              Z_in  = Z_in  + Z(el_in,dof)*phie(dof,gp_in,et)
              Qx_in = Qx_in + Qx(el_in,dof)*phie(dof,gp_in,et)
              Qy_in = Qy_in + Qy(el_in,dof)*phie(dof,gp_in,et)
            ENDDO            
              
            Z_ex = Z_in
            Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
            Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

            CALL numerical_flux(Qx_in,Qy_in,Z_in,Qx_ex,Qy_ex,Z_ex,hb,nx,ny,sp,Qxhat,Qyhat,Zhat)
           
            Zhat  = detJe(ged,pt)*Zhat              
            Qxhat = detJe(ged,pt)*Qxhat
            Qyhat = detJe(ged,pt)*Qyhat
              
            DO l = 1,ndof(et)
              rhsZ(el_in,l)  = rhsZ(el_in,l)  - Zhat*phie_int(l,gp_in,et)
              rhsQx(el_in,l) = rhsQx(el_in,l) - Qxhat*phie_int(l,gp_in,et)
              rhsQy(el_in,l) = rhsQy(el_in,l) - Qyhat*phie_int(l,gp_in,et)                   
            ENDDO  
              
          ENDDO
        ENDDO
!$OMP end do
      ENDIF     


      END SUBROUTINE boundary_edge_land