       SUBROUTINE interior_edge_calculations(et,sel,eel,sed,eed,ndof,nqpte)

       USE globals, ONLY: Hqpt,Qxqpt,Qyqpt, &
                          H,Qx,Qy, &
                          phie, &
                          g,pt5g,recipHa,xmom,ymom,xymom, &
                          const,inx,iny, &
                          Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Hhatv,Qxhatv,Qyhatv, &
                          inx,iny,detJe_in,detJe_ex, &
                          Hfe,Hfi,Qxfi,Qxfe,Qyfi,Qyfe                                                 

       IMPLICIT NONE
       
       INTEGER :: pt,el,dof,ed,et
       INTEGER :: nqpte,ndof
       INTEGER :: sel,eel
       INTEGER :: sed,eed
       
ed_points: DO pt = 1,3*nqpte

!!DIR$ VECTOR ALIGNED               
              DO el = sel,eel
                Hqpt(el,pt) = H(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof     
!!DIR$ VECTOR ALIGNED
                DO el = sel,eel   ! Compute solutions at edge quadrature points                
                  Hqpt(el,pt) = Hqpt(el,pt) + H(el,dof)*phie(dof,pt,et)            
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,et)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,et)            
                ENDDO

              ENDDO ed_basis

!!DIR$ VECTOR ALIGNED
              DO el = sel,eel  ! Compute momentum terms       
                recipHa(el) = 1d0/Hqpt(el,pt)
                
                xmom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xymom(el,pt) = Qxqpt(el,pt)*Qyqpt(el,pt)*recipHa(el)
              ENDDO

          ENDDO ed_points
                  
        
ed_points2: DO pt = 1,nqpte ! Compute numerical fluxes for all edges
              
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO
                                        
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed            
                Hhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
              ENDDO
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed        
                Qxhatv(ed) = .5d0*(inx(ed,pt)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))
              ENDDO
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                Qyhatv(ed) = .5d0*(inx(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed,pt)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))
              ENDDO
!DIR$ IVDEP              
              DO ed = sed,eed 
                Hfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  detJe_in(ed,pt)*Hhatv(ed)
              ENDDO          
!DIR$ IVDEP              
              DO ed = sed,eed                                
                Qxfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP              
              DO ed = sed,eed                                    
                Qyfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qyhatv(ed)
              ENDDO               

        ENDDO ed_points2       


       RETURN 
       END SUBROUTINE
       
       
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
       
       
       
       SUBROUTINE edge_integration(et,sel,eel,ndof,nqpte)
       
       USE globals, ONLY: rhsH,rhsQx,rhsQy, &
                          Hflux,Qxflux,Qyflux, &
                          phie_int
       
       IMPLICIT NONE
       
       INTEGER pt,l,el
       INTEGER :: et
       INTEGER :: sel,eel
       INTEGER :: ndof,nqpte
       
       DO pt = 1,3*nqpte
         DO l = 1,ndof
!!DIR$ VECTOR ALIGNED
           DO el = sel,eel
             rhsH(el,l) = rhsH(el,l) - Hflux(el,pt)*phie_int(l,pt,et)
             rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt,et)
             rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt,et)                   
           ENDDO
         ENDDO                                    
       ENDDO       
       
       RETURN
       END SUBROUTINE