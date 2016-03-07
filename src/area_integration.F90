      SUBROUTINE area_integration(et,sel,eel,ndof,nqpta)
      
      USE globals, ONLY: Zqpt,Hqpt,Qxqpt,Qyqpt, &
                         Z,Qx,Qy, &
                         phia,phia_int,dpdx,dpdy, &
                         recipHa,xmom,ymom,xymom, &
                         tau,src_x,src_y, &
                         dhbdx,dhbdy, &
                         rhsZ,rhsQx,rhsQy,&
                         g,pt5g,hbqpta
                         
      USE read_dginp, ONLY: cf                         

      IMPLICIT NONE

      INTEGER :: nqpta,ndof,et
      INTEGER :: sel,eel

      INTEGER :: pt,dof,el,l
      INTEGER :: ind


a_points: DO pt = 1,nqpta

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel  ! First basis function is 1
              Zqpt(el,1)  = Z(el,1)
              Qxqpt(el,1) = Qx(el,1)
              Qyqpt(el,1) = Qy(el,1)
            ENDDO

   a_basis: DO dof = 2,ndof
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
              DO el = sel,eel   ! Evaluate solution at area quadrature point
                Zqpt(el,1)  = Zqpt(el,1)  + Z(el,dof)*phia(dof,pt,et)
                Qxqpt(el,1) = Qxqpt(el,1) + Qx(el,dof)*phia(dof,pt,et) 
                Qyqpt(el,1) = Qyqpt(el,1) + Qy(el,dof)*phia(dof,pt,et)
              ENDDO

            ENDDO a_basis

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Compute momentum terms
              Hqpt(el,1) = Zqpt(el,1) + hbqpta(el,pt)
              recipHa(el) = 1d0/Hqpt(el,1)

              xmom(el,1) = pt5g*(Hqpt(el,1)*Hqpt(el,1)-hbqpta(el,pt)*hbqpta(el,pt)) + Qxqpt(el,1)*Qxqpt(el,1)*recipHa(el)
              ymom(el,1) = pt5g*(Hqpt(el,1)*Hqpt(el,1)-hbqpta(el,pt)*hbqpta(el,pt)) + Qyqpt(el,1)*Qyqpt(el,1)*recipHa(el) 
              xymom(el,1) = Qxqpt(el,1)*Qyqpt(el,1)*recipHa(el)
            ENDDO 

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el =  sel,eel   ! Compute source terms
              tau(el) = cf*sqrt((Qxqpt(el,1)*recipHa(el))**2 + (Qyqpt(el,1)*recipHa(el))**2)*recipHa(el)
              src_x(el) = g*Zqpt(el,1)*dhbdx(el,pt) - tau(el)*Qxqpt(el,1) 
              src_y(el) = g*Zqpt(el,1)*dhbdy(el,pt) - tau(el)*Qyqpt(el,1)
            ENDDO

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Derivatives are 0 for first dof
              rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(el,1,pt)
              rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(el,1,pt)
            ENDDO

      test: DO l = 2,ndof 
!!DIR$ VECTOR ALIGNED          
!DIR$ SIMD
              DO el = sel,eel
                rhsZ(el,l)  = rhsZ(el,l)  + Qxqpt(el,1)*dpdx(el,l,pt) + Qyqpt(el,1)*dpdy(el,l,pt)

                rhsQx(el,l) = rhsQx(el,l) + xmom(el,1)*dpdx(el,l,pt)  + xymom(el,1)*dpdy(el,l,pt) + src_x(el)*phia_int(el,l,pt)           
                rhsQy(el,l) = rhsQy(el,l) + xymom(el,1)*dpdx(el,l,pt) + ymom(el,1)*dpdy(el,l,pt)  + src_y(el)*phia_int(el,l,pt)
              ENDDO

            ENDDO test 

          ENDDO a_points
          
      RETURN 
      END SUBROUTINE