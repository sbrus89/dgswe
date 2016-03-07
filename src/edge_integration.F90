       SUBROUTINE interior_edge_eval(et,sel,eel,ndof,tnqpte)

       USE globals, ONLY: Zqpt,Hqpt,Qxqpt,Qyqpt, &
                          Z,Qx,Qy, &
                          phie,hbqpte, &
                          g,pt5g,recipHa,xmom,ymom,xymom
                                               

       IMPLICIT NONE
       
       INTEGER :: pt,el,dof,ed,et
       INTEGER :: tnqpte,ndof
       INTEGER :: sel,eel
       
ed_points: DO pt = 1,tnqpte

!!DIR$ VECTOR ALIGNED     
!DIR$ SIMD
              DO el = sel,eel
                Zqpt(el,pt)  = Z(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof     
!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
                DO el = sel,eel   ! Compute solutions at edge quadrature points                
                  Zqpt(el,pt)  = Zqpt(el,pt)  + Z(el,dof)*phie(dof,pt,et)            
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,et)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,et)            
                ENDDO

              ENDDO ed_basis

!!DIR$ VECTOR ALIGNED
!DIR$ SIMD
              DO el = sel,eel  ! Compute momentum terms   
                Hqpt(el,pt) = Zqpt(el,pt) + hbqpte(el,pt)    
                recipHa(el) = 1d0/Hqpt(el,pt)
                
                xmom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xymom(el,pt) = Qxqpt(el,pt)*Qyqpt(el,pt)*recipHa(el)
              ENDDO

          ENDDO ed_points

      END SUBROUTINE interior_edge_eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE interior_edge_nflux(sed,eed,nqpte)

      USE globals, ONLY:  const,inx,iny,icfac,g, &
                          Zi,Ze,Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Zhatv,Qxhatv,Qyhatv, &
                          detJe_in,detJe_ex, &
                          Zfe,Zfi,Qxfi,Qxfe,Qyfi,Qyfe  
   
      IMPLICIT NONE

      INTEGER :: sed,eed,nqpte
      INTEGER :: pt,ed
                  
        
ed_points: DO pt = 1,nqpte ! Compute numerical fluxes for all edges
              
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr*icfac(ed,pt)), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr*icfac(ed,pt)))
              ENDDO
                                        
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed                  
                Zhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(Ze(ed,pt)%ptr - Zi(ed,pt)%ptr))           
              ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed     
                Qxhatv(ed) = .5d0*(inx(ed,pt)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))                                     
              ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
              DO ed = sed,eed
                Qyhatv(ed) = .5d0*(inx(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed,pt)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))                                     
              ENDDO
!DIR$ IVDEP              
              DO ed = sed,eed
                Zfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Zhatv(ed)
                Zfi(ed,pt)%ptr =  detJe_in(ed,pt)*Zhatv(ed)                
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

        ENDDO ed_points  


       RETURN 
       END SUBROUTINE interior_edge_nflux
       
       
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
      SUBROUTINE recieve_edge_nflux(nred,nqpte)

      USE globals, ONLY: const,g,pt5g,recipHa, &
                         Zhatv,Qxhatv,Qyhatv
         
      USE messenger2, ONLY: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre, &
                            xmri,xmre,ymri,ymre,xymri,xymre, &
                            detJe_recv,rcfac,hbr,rnx,rny, &
                            Zfri,Qxfri,Qyfri  
   
      IMPLICIT NONE

      INTEGER :: nred,nqpte
      INTEGER :: pt,ed
                  
        
      DO pt = 1,nqpte
      
        DO ed = 1,nred
          Hre(ed) = Zre(ed,pt)%ptr + hbr(ed,pt)
        ENDDO
      
!!DIR$ VECTOR ALIGNED      
        DO ed = 1,nred
          const(ed) = max(abs(Qxri(ed,pt)%ptr*rnx(ed,pt) + Qyri(ed,pt)%ptr*rny(ed,pt))/Hri(ed,pt)%ptr + sqrt(g*Hri(ed,pt)%ptr*rcfac(ed,pt)), &
                          abs(Qxre(ed,pt)%ptr*rnx(ed,pt) + Qyre(ed,pt)%ptr*rny(ed,pt))/Hre(ed)        + sqrt(g*Hre(ed)*rcfac(ed,pt)))          
        ENDDO
        

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
        DO ed = 1,nred      
          Zhatv(ed) = .5d0*(rnx(ed,pt)*(Qxri(ed,pt)%ptr + Qxre(ed,pt)%ptr) + rny(ed,pt)*(Qyri(ed,pt)%ptr + Qyre(ed,pt)%ptr) &
                                        - const(ed)*(Zre(ed,pt)%ptr - Zri(ed,pt)%ptr))
        ENDDO

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED        
        DO ed = 1,nred
          recipHa(ed) = 1d0/Hre(ed)
          
          xmre(ed) = pt5g*(Hre(ed)*Hre(ed) - hbr(ed,pt)*hbr(ed,pt)) + Qxre(ed,pt)%ptr*Qxre(ed,pt)%ptr*recipHa(ed)
          ymre(ed) = pt5g*(Hre(ed)*Hre(ed) - hbr(ed,pt)*hbr(ed,pt)) + Qyre(ed,pt)%ptr*Qyre(ed,pt)%ptr*recipHa(ed)
          xymre(ed) = Qxre(ed,pt)%ptr*Qyre(ed,pt)%ptr*recipHa(ed)
        ENDDO
 
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED 
        DO ed = 1,nred
          Qxhatv(ed) = .5d0*(rnx(ed,pt)*(xmri(ed,pt)%ptr + xmre(ed)) + rny(ed,pt)*(xymri(ed,pt)%ptr + xymre(ed))  &
                                        - const(ed)*(Qxre(ed,pt)%ptr - Qxri(ed,pt)%ptr))
        ENDDO
  
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED  
        DO ed = 1,nred
          Qyhatv(ed) = .5d0*(rnx(ed,pt)*(xymri(ed,pt)%ptr + xymre(ed)) + rny(ed,pt)*(ymri(ed,pt)%ptr + ymre(ed))  &
                                        - const(ed)*(Qyre(ed,pt)%ptr - Qyri(ed,pt)%ptr))
        ENDDO

        DO ed = 1,nred                                       
          Zfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Zhatv(ed)          
        ENDDO   
        DO ed = 1,nred                                         
          Qxfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qxhatv(ed)        
        ENDDO        
        DO ed = 1,nred 
          Qyfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qyhatv(ed)        
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE recieve_edge_nflux
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
       
       
       SUBROUTINE edge_integration(et,sel,eel,ndof,tnqpte)
       
       USE globals, ONLY: rhsZ,rhsQx,rhsQy, &
                          Zflux,Qxflux,Qyflux, &
                          phie_int
       
       IMPLICIT NONE
       
       INTEGER pt,l,el
       INTEGER :: et
       INTEGER :: sel,eel
       INTEGER :: ndof,tnqpte
       
       DO pt = 1,tnqpte
         DO l = 1,ndof
!!DIR$ VECTOR ALIGNED
           DO el = sel,eel
             rhsZ(el,l) = rhsZ(el,l) - Zflux(el,pt)*phie_int(l,pt,et)
             rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt,et)
             rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt,et)                   
           ENDDO
         ENDDO                                    
       ENDDO       
       
       RETURN
       END SUBROUTINE