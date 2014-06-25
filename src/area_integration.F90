      SUBROUTINE area_integration(et,sel,eel,ndof,nqpta)
      
      USE globals, ONLY: Hqpt,Qxqpt,Qyqpt, &
                         H,Qx,Qy, &
                         phia,phia_int,dpdx,dpdy, &
                         recipHa,xmom,ymom,xymom, &
                         tau,src_x,src_y, &
                         dhbdx,dhbdy, &
                         rhsH,rhsQx,rhsQy,&
                         g,pt5g,cf

      IMPLICIT NONE

      INTEGER :: pt,dof,el,l
      INTEGER :: ind
      INTEGER :: nqpta,ndof,et
      INTEGER :: sel,eel



a_points: DO pt = 1,nqpta

!!DIR$ VECTOR ALIGNED
            DO el = sel,eel  ! First basis function is 1
              Hqpt(el,1)  = H(el,1)
              Qxqpt(el,1) = Qx(el,1)
              Qyqpt(el,1) = Qy(el,1)
            ENDDO

   a_basis: DO dof = 2,ndof
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
              DO el = sel,eel   ! Evaluate solution at area quadrature point
                Hqpt(el,1)  = Hqpt(el,1)  + H(el,dof)*phia(dof,pt,et)
                Qxqpt(el,1) = Qxqpt(el,1) + Qx(el,dof)*phia(dof,pt,et) 
                Qyqpt(el,1) = Qyqpt(el,1) + Qy(el,dof)*phia(dof,pt,et)
              ENDDO

            ENDDO a_basis

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Compute momentum terms
              recipHa(el) = 1d0/Hqpt(el,1)

              xmom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qxqpt(el,1)*Qxqpt(el,1)*recipHa(el)
              ymom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qyqpt(el,1)*Qyqpt(el,1)*recipHa(el) 
              xymom(el,1) = Qxqpt(el,1)*Qyqpt(el,1)*recipHa(el)
            ENDDO 

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el =  sel,eel   ! Compute source terms
              tau(el) = cf*sqrt((Qxqpt(el,1)*recipHa(el))**2 + (Qyqpt(el,1)*recipHa(el))**2)*recipHa(el)
              src_x(el) = g*Hqpt(el,1)*dhbdx(el,pt) - tau(el)*Qxqpt(el,1) 
              src_y(el) = g*Hqpt(el,1)*dhbdy(el,pt) - tau(el)*Qyqpt(el,1)
            ENDDO

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
            DO el = sel,eel   ! Derivatives are 0 for first dof
              rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(el,pt)
              rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(el,pt)
            ENDDO

      test: DO l = 2,ndof 
              ind = (l-1)*nqpta+pt
!!DIR$ VECTOR ALIGNED          
!DIR$ SIMD
              DO el = sel,eel
                rhsH(el,l)  = rhsH(el,l)  + Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind)

                rhsQx(el,l) = rhsQx(el,l) + xmom(el,1)*dpdx(el,ind)  + xymom(el,1)*dpdy(el,ind) + src_x(el)*phia_int(el,ind)           
                rhsQy(el,l) = rhsQy(el,l) + xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind)  + src_y(el)*phia_int(el,ind)
              ENDDO

            ENDDO test 

          ENDDO a_points
          
      RETURN 
      END SUBROUTINE