
      SUBROUTINE interior_edge_nflux(spt,ept)

      USE globals, ONLY:  const,inx,iny,icfac,g, &
                          Zi,Ze,Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Zhatv,Qxhatv,Qyhatv, &
                          detJe_in,detJe_ex, &
                          Exxi,Exxe,Eyyi,Eyye,Exyi,Exye,Eyxi,Eyxe, &
                          esl_tol
                          
      USE read_dginp, ONLY: esl
   
      IMPLICIT NONE

      INTEGER :: spt,ept
      INTEGER :: pt
                  
        
              
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
      DO pt = spt,ept
                
        const(pt) = max(abs(Qxi(pt)%ptr*inx(pt) + Qyi(pt)%ptr*iny(pt))/Hi(pt)%ptr &
                                                      + sqrt(g*Hi(pt)%ptr*icfac(pt)), &
                        abs(Qxe(pt)%ptr*inx(pt) + Qye(pt)%ptr*iny(pt))/He(pt)%ptr &
                                                      + sqrt(g*He(pt)%ptr*icfac(pt)))
      ENDDO
                                        
!DIR$ SIMD              
!!DIR$ VECTOR ALIGNED
      DO pt = spt,ept                
        Zhatv(pt) = .5d0*(inx(pt)*(Qxi(pt)%ptr + Qxe(pt)%ptr) + iny(pt)*(Qyi(pt)%ptr + Qye(pt)%ptr) &
                              - const(pt)*(Ze(pt)%ptr - Zi(pt)%ptr))           
      ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
      DO pt = spt,ept     
        Qxhatv(pt) = .5d0*(inx(pt)*(xmi(pt)%ptr + xme(pt)%ptr) + iny(pt)*(xymi(pt)%ptr + xyme(pt)%ptr)  &
                               - const(pt)*(Qxe(pt)%ptr - Qxi(pt)%ptr))                                     
      ENDDO
              
!DIR$ SIMD
!!DIR$ VECTOR ALIGNED
      DO pt = spt,ept
        Qyhatv(pt) = .5d0*(inx(pt)*(xymi(pt)%ptr + xyme(pt)%ptr) + iny(pt)*(ymi(pt)%ptr + yme(pt)%ptr)  &
                               - const(pt)*(Qye(pt)%ptr - Qyi(pt)%ptr))                                     
      ENDDO
              
              
              
              
      ! Eddy viscosity terms
      IF (esl > esl_tol) THEN
        DO pt = spt,ept
          Qxhatv(pt) = Qxhatv(pt) - inx(pt)*.5d0*(Exxe(pt)%ptr+Exxi(pt)%ptr) &
                                  - iny(pt)*.5d0*(Exye(pt)%ptr+Exyi(pt)%ptr)               
        ENDDO
              
        DO pt = spt,ept
          Qyhatv(pt) = Qyhatv(pt) - inx(pt)*.5d0*(Exye(pt)%ptr+Exyi(pt)%ptr) &
                                  - iny(pt)*.5d0*(Eyye(pt)%ptr+Eyyi(pt)%ptr)                
        ENDDO
      ENDIF
              
              
              
              
!DIR$ IVDEP              
      DO pt = spt,ept
        Ze(pt)%ptr = -detJe_ex(pt)*Zhatv(pt)
        Zi(pt)%ptr =  detJe_in(pt)*Zhatv(pt)                
      ENDDO          
!DIR$ IVDEP              
      DO pt = spt,ept                
        Qxe(pt)%ptr = -detJe_ex(pt)*Qxhatv(pt)
        Qxi(pt)%ptr =  detJe_in(pt)*Qxhatv(pt)
      ENDDO   
!DIR$ IVDEP              
      DO pt = spt,ept                      
        Qye(pt)%ptr = -detJe_ex(pt)*Qyhatv(pt)
        Qyi(pt)%ptr =  detJe_in(pt)*Qyhatv(pt)
      ENDDO               


      RETURN 
      END SUBROUTINE interior_edge_nflux
       
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE interior_edge_nflux_ldg(spt,ept)  
      
      USE globals, ONLY: Exxi,Exxe,Eyyi,Eyye,Exyi,Exye,Eyxi,Eyxe, &
                         Qxi,Qxe,Qyi,Qye, &
                         inx,iny, &
                         detJe_in,detJe_ex
      
      IMPLICIT NONE
      
      INTEGER :: spt,ept
      INTEGER :: nqpte
      
      INTEGER :: pt,ed
                  
      DO pt = spt,ept

          Exxe(pt)%ptr = -detJe_ex(pt)*inx(pt)*.5d0*(Qxi(pt)%ptr+Qxe(pt)%ptr)
          Exxi(pt)%ptr =  detJe_in(pt)*inx(pt)*.5d0*(Qxe(pt)%ptr+Qxi(pt)%ptr)

          Eyye(pt)%ptr = -detJe_ex(pt)*iny(pt)*.5d0*(Qyi(pt)%ptr+Qye(pt)%ptr)           
          Eyyi(pt)%ptr =  detJe_in(pt)*iny(pt)*.5d0*(Qye(pt)%ptr+Qyi(pt)%ptr)
          
          Exye(pt)%ptr = -detJe_ex(pt)*.5d0*(iny(pt)*(Qxi(pt)%ptr+Qxe(pt)%ptr) &
                                            +inx(pt)*(Qyi(pt)%ptr+Qye(pt)%ptr))
          Exyi(pt)%ptr =  detJe_in(pt)*.5d0*(iny(pt)*(Qxe(pt)%ptr+Qxi(pt)%ptr) &
                                            +inx(pt)*(Qye(pt)%ptr+Qyi(pt)%ptr))          
      ENDDO
      
      RETURN
      END SUBROUTINE interior_edge_nflux_ldg
                            
                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
      SUBROUTINE recieve_edge_nflux(nred,nqpte)

      USE globals, ONLY: const,g,pt5g,recipHa, &
                         Zhatv,Qxhatv,Qyhatv, &
                         esl_tol
         
      USE messenger2, ONLY: Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre, &
                            xmri,xmre,ymri,ymre,xymri,xymre, &
                            detJe_recv,rcfac,hbr,rnx,rny, &
                            Exxri,Exxre,Eyyri,Eyyre,Exyri,Exyre,Eyxri,Eyxre
                            
      USE read_dginp, ONLY: esl
   
      IMPLICIT NONE

      INTEGER :: nred,nqpte
      INTEGER :: pt,ed
                  
        
      DO pt = 1,nqpte
      
        DO ed = 1,nred
          Hre(ed) = Zre(ed,pt)%ptr + hbr(ed,pt)
        ENDDO
      
!!DIR$ VECTOR ALIGNED      
        DO ed = 1,nred
          const(ed) = max(abs(Qxri(ed,pt)%ptr*rnx(ed,pt) + Qyri(ed,pt)%ptr*rny(ed,pt))/Hri(ed,pt)%ptr  &
                                                         + sqrt(g*Hri(ed,pt)%ptr*rcfac(ed,pt)), &
                          abs(Qxre(ed,pt)%ptr*rnx(ed,pt) + Qyre(ed,pt)%ptr*rny(ed,pt))/Hre(ed) &
                                                         + sqrt(g*Hre(ed)*rcfac(ed,pt)))          
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
        
        
        ! Eddy viscosity terms
        IF (esl > esl_tol) THEN
          DO ed = 1,nred
            Qxhatv(ed) = Qxhatv(ed) - rnx(ed,pt)*.5d0*(Exxri(ed,pt)%ptr+Exxre(ed,pt)%ptr)  &
                                    - rny(ed,pt)*.5d0*(Exyri(ed,pt)%ptr+Exyre(ed,pt)%ptr)
          ENDDO
        
          DO ed = 1,nred
            Qyhatv(ed) = Qyhatv(ed) - rnx(ed,pt)*.5d0*(Exyri(ed,pt)%ptr+Exyre(ed,pt)%ptr)  &
                                    - rny(ed,pt)*.5d0*(Eyyri(ed,pt)%ptr+Eyyre(ed,pt)%ptr)          
          ENDDO
        ENDIF
        
        
        DO ed = 1,nred                                       
          Zri(ed,pt)%ptr =  detJe_recv(ed,pt)*Zhatv(ed)          
        ENDDO   
        DO ed = 1,nred                                            
          Qxri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qxhatv(ed)
        ENDDO        
        DO ed = 1,nred 
          Qyri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qyhatv(ed)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE recieve_edge_nflux
       
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE recieve_edge_nflux_ldg(nred,nqpte)

      USE messenger2, ONLY: Qxre,Qyre,Qxri,Qyri, &
                            Exxri,Eyyri,Exyri,Eyxri, &
                            detJe_recv,rnx,rny
      
      IMPLICIT NONE
      
      INTEGER :: nred,nqpte
      
      INTEGER :: pt,ed
      
      DO pt = 1,nqpte
        DO ed = 1,nred

          Exxri(ed,pt)%ptr = detJe_recv(ed,pt)*rnx(ed,pt)*.5d0*(Qxre(ed,pt)%ptr+Qxri(ed,pt)%ptr)
          Eyyri(ed,pt)%ptr = detJe_recv(ed,pt)*rny(ed,pt)*.5d0*(Qyre(ed,pt)%ptr+Qyri(ed,pt)%ptr)
          Exyri(ed,pt)%ptr = detJe_recv(ed,pt)*.5d0*(rny(ed,pt)*(Qxre(ed,pt)%ptr+Qxri(ed,pt)%ptr) &
                                                    +rnx(ed,pt)*(Qyre(ed,pt)%ptr+Qyri(ed,pt)%ptr))
          
        ENDDO
      ENDDO
                  
      RETURN
      END SUBROUTINE recieve_edge_nflux_ldg

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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