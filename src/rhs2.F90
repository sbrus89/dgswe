      SUBROUTINE rhs2()

      USE globals,  ONLY: pres,l,ind,el,ne,pt,led,ed,dof,ndof,nel_type,el_type, &
                          blk,elblk,nfblk,rnfblk,nrblk,npartet,nverts, &
                          g,pt5g,tstage,ramp,arg,bfr, & 
                          H,Z,Qx,Qy, &
                          rhsH,rhsZ,rhsQx,rhsQy, &
                          Hqpt,Zqpt,Qxqpt,Qyqpt, &
                          hbqpta,hbqpte,hbqpted, &
                          nqpta,wpta,phia,phia_int,dpdx,dpdy,nqpte,wpte,phie,phie_int,mmi, &
                          press,recipH,xmom,ymom,xymom,tau,u,v,src_x,src_y,dhbdx,dhbdy, &
                          el_in,el_ex,led_in,led_ex,gp_in,gp_ex,check_iedge,check_gedge, &
                          nx,ny,nx2,ny2,nxny,tx,ty, &
                          H_in,H_ex,Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex,Qn,Qt, &
                          xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                          Hhat,Zhat,Qxhat,Qyhat, &
                          Hflux,Zflux,Qxflux,Qyflux, &
                          nied,iedn,nobed,obedn,nnfbed,nfbedn,nfbed,fbedn, &
                          nobfr,obfreq,obper,obnfact,obamp_qpt,obph_qpt,obeq,obdepth_qpt, &
                          nfbfr,fbfreq,fbper,fbnfact,fbamp_qpt,fbph_qpt,fbeq, &
                          ged2led,ged2el,gel2ael,ged, &
                          pressa,recipHa, & 
                          Hi,He,Zi,Ze,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Hfi,Hfe,Zfi,Zfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                          const,inx,iny,icfac,detJe_in,detJe_ex,detJe, &
                          nx_pt,ny_pt,Spe, &
                          Hhatv,Zhatv,Qxhatv,Qyhatv, &
                          MirhsH,MirhsZ,MirhsQx,MirhsQy
                          
      USE read_dginp, ONLY: npart,cf                          
                          
#ifdef CMPI                          
      USE messenger2, ONLY: myrank,message_recieve,message_send, &
                            nred,nproc_sr,match_edge, &
                            solreq,solreq_send,solreq_recv,ierr, &
                            Zri,Zre,Hri,Hre,Qxri,Qxre,Qyri,Qyre, &
                            xmri,ymri,xymri,xmre,ymre,xymre, &
                            Zfri,Hfri,Qxfri,Qyfri, &
                            rnx,rny,rcfac,detJe_recv,hbr         
                            
      USE mpi                            
#endif       
!$    USE omp_lib                          
                          
                   
      IMPLICIT NONE
      INTEGER :: i,j,et,m,ete
      REAL(pres) :: sp,hb

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

!!DIR$ VECTOR ALIGNED
!       rhsH(:,:) = 0d0
      rhsZ(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQy(:,:) = 0d0   
      
!!DIR$ VECTOR ALIGNED
!       MirhsH(:,:) = 0d0
      MirhsZ(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      MirhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      MirhsQy(:,:) = 0d0         

      
!$OMP parallel default(none)  &
!$OMP             private(ind,blk,pt,dof,el,et,ete,l,m,ed,ged,led_in,gp_in,el_in, &
!$OMP                     nx,ny,tx,ty,nx2,ny2,nxny,sp,hb,H_in,H_ex,Z_in,Z_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!$OMP                     press,recipH,xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!$OMP                     Hhat,Zhat,Qxhat,Qyhat, &
!$OMP                     Qn,Qt,arg,bfr) &
!$OMP             shared(nqpta,nqpte,ndof,nverts,H,Z,Qx,Qy,Hqpt,Qxqpt,Qyqpt, & 
!$OMP                    phia,phie,recipHa,xmom,ymom,xymom,tau,src_x,src_y, &
!$OMP                    dhbdx,dhbdy,cf,phia_int,phie_int,dpdx,dpdy,mmi, &
!$OMP                    MirhsH,MirhsZ,MirhsQx,MirhsQy,rhsH,rhsQx,rhsQy, &
!$OMP                    const,Qxi,Qxe,Qyi,Qye,Hi,He,Zi,Ze,xmi,xme,ymi,yme,xymi,xyme, &
!$OMP                    Hhatv,Zhatv,Qxhatv,Qyhatv,Hfi,Hfe,Zfi,Zfe,Qxfi,Qxfe,Qyfi,Qyfe, &
!$OMP                    Hflux,Zflux,Qxflux,Qyflux,inx,iny,icfac,detJe_in,detJe_ex, &
!$OMP                    nnfbed,nfbedn,nfbed,fbedn,nobed,obedn, &
!$OMP                    nfbfr,fbfreq,fbper,fbeq,fbamp_qpt,fbnfact,fbph_qpt, &
!$OMP                    nobfr,obfreq,obper,obeq,obamp_qpt,obnfact,obph_qpt,obdepth_qpt, &
!$OMP                    ramp,tstage, &
!$OMP                    detJe,nx_pt,ny_pt,Spe, &
!$OMP                    ged2led,gel2ael,ged2el,el_type, &
!$OMP                    npart,nrblk,elblk,nfblk,rnfblk,npartet)


!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

!$OMP do
   

      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
!             CALL area_integration(et,elblk(1,blk,et),elblk(2,blk,et),ndof(et),nqpta(et))
            
a_points:   DO pt = 1,nqpta(et)

!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)  ! First basis function is 1
!                 Hqpt(el,1)  = H(el,1)
                Zqpt(el,1)  = Z(el,1)
                Qxqpt(el,1) = Qx(el,1)
                Qyqpt(el,1) = Qy(el,1)
              ENDDO

              DO dof = 2,ndof(et)
!!DIR$ VECTOR ALIGNED
! ! DIR$ SIMD
                DO el = elblk(1,blk,et),elblk(2,blk,et)   ! Evaluate solution at area quadrature point
!                   Hqpt(el,1)  = Hqpt(el,1)  + H(el,dof)*phia(dof,pt,et)
                  Zqpt(el,1)  = Zqpt(el,1)  + Z(el,dof)*phia(dof,pt,et)
                  Qxqpt(el,1) = Qxqpt(el,1) + Qx(el,dof)*phia(dof,pt,et) 
                  Qyqpt(el,1) = Qyqpt(el,1) + Qy(el,dof)*phia(dof,pt,et)
                ENDDO
              ENDDO 

!!DIR$ VECTOR ALIGNED
! !DIR$ SIMD
              DO el = elblk(1,blk,et),elblk(2,blk,et)   ! Compute momentum terms
                Hqpt(el,1) = Zqpt(el,1) + hbqpta(el,pt)
                
                recipHa(el) = 1d0/Hqpt(el,1)

!                 xmom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qxqpt(el,1)*Qxqpt(el,1)*recipHa(el)
!                 ymom(el,1) = pt5g*Hqpt(el,1)*Hqpt(el,1) + Qyqpt(el,1)*Qyqpt(el,1)*recipHa(el) 

                xmom(el,1) = pt5g*(Hqpt(el,1)*Hqpt(el,1)-hbqpta(el,pt)*hbqpta(el,pt)) + Qxqpt(el,1)*Qxqpt(el,1)*recipHa(el)
                ymom(el,1) = pt5g*(Hqpt(el,1)*Hqpt(el,1)-hbqpta(el,pt)*hbqpta(el,pt)) + Qyqpt(el,1)*Qyqpt(el,1)*recipHa(el) 
                xymom(el,1) = Qxqpt(el,1)*Qyqpt(el,1)*recipHa(el)
              ENDDO 

!!DIR$ VECTOR ALIGNED
! !DIR$ SIMD
              DO el =  elblk(1,blk,et),elblk(2,blk,et)   ! Compute source terms
                tau(el) = cf*sqrt((Qxqpt(el,1)*recipHa(el))**2 + (Qyqpt(el,1)*recipHa(el))**2)*recipHa(el)
!                 src_x(el) = g*Hqpt(el,1)*dhbdx(el,pt) - tau(el)*Qxqpt(el,1) 
!                 src_y(el) = g*Hqpt(el,1)*dhbdy(el,pt) - tau(el)*Qyqpt(el,1)

                src_x(el) = g*Zqpt(el,1)*dhbdx(el,pt) - tau(el)*Qxqpt(el,1) 
                src_y(el) = g*Zqpt(el,1)*dhbdy(el,pt) - tau(el)*Qyqpt(el,1)
              ENDDO

!!DIR$ VECTOR ALIGNED
! !DIR$ SIMD
              DO el = elblk(1,blk,et),elblk(2,blk,et)   ! Derivatives are 0 for first dof
                rhsQx(el,1) = rhsQx(el,1) + src_x(el)*phia_int(el,pt)
                rhsQy(el,1) = rhsQy(el,1) + src_y(el)*phia_int(el,pt)
              ENDDO

             DO l = 2,ndof(et) 
                ind = (l-1)*nqpta(et)+pt
!!DIR$ VECTOR ALIGNED          
! !DIR$ SIMD
                DO el = elblk(1,blk,et),elblk(2,blk,et)
!                   rhsH(el,l)  = rhsH(el,l)  + Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind)

                  rhsZ(el,l)  = rhsZ(el,l)  + Qxqpt(el,1)*dpdx(el,ind) + Qyqpt(el,1)*dpdy(el,ind)

                  rhsQx(el,l) = rhsQx(el,l) + xmom(el,1)*dpdx(el,ind)  + xymom(el,1)*dpdy(el,ind) + src_x(el)*phia_int(el,ind)           
                  rhsQy(el,l) = rhsQy(el,l) + xymom(el,1)*dpdx(el,ind) + ymom(el,1)*dpdy(el,ind)  + src_y(el)*phia_int(el,ind)
                ENDDO
              ENDDO  

            ENDDO a_points 
            
          ENDIF
        ENDDO
      ENDDO      
      
      
!$OMP end do 
      
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Edge Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

#ifdef CMPI 
!       post a non-blocking recieve from all processes
      CALL message_recieve()
#endif    
      

!$OMP do
!       DO blk = 1,npart
!         DO et = 1,nel_type       
!           IF (npartet(et,blk) > 0) THEN
!             CALL interior_edge_calculations(et,elblk(1,blk,et),elblk(2,blk,et),nfblk(1,blk),nfblk(2,blk),ndof(et),nqpta(et))
!           ENDIF
!         ENDDO
!       ENDDO
      
       DO blk = 1,npart 
        DO et = 1,nel_type   
        
          IF (et == 3) THEN 
            ete = 1
          ELSE IF (et == 4) THEN
            ete = 2
          ELSE   
            ete = et
          ENDIF
          
          IF (npartet(et,blk) > 0) THEN       
ed_points: DO pt = 1,nverts(et)*nqpte(ete)

!!DIR$ VECTOR ALIGNED               
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 Hqpt(el,pt) = H(el,1)
                Zqpt(el,pt)  = Z(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof(et)
!!DIR$ VECTOR ALIGNED
                DO el = elblk(1,blk,et),elblk(2,blk,et)   ! Compute solutions at edge quadrature points                
!                   Hqpt(el,pt) = Hqpt(el,pt) + H(el,dof)*phie(dof,pt,ete)   
                  Zqpt(el,pt)  = Zqpt(el,pt)  + Z(el,dof)*phie(dof,pt,ete)                  
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,ete)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,ete)            
                ENDDO

              ENDDO ed_basis

!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)  ! Compute momentum terms  
                Hqpt(el,pt) = Zqpt(el,pt) + hbqpte(el,pt)
                recipHa(el) = 1d0/Hqpt(el,pt)
                
!                 xmom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
!                 ymom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xmom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*(Hqpt(el,pt)*Hqpt(el,pt) - hbqpte(el,pt)*hbqpte(el,pt))+ Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
                xymom(el,pt) = Qxqpt(el,pt)*Qyqpt(el,pt)*recipHa(el)
              ENDDO

          ENDDO ed_points
         ENDIF
        ENDDO
        
#ifdef CMPI    
       ! Post an non-blocking send to all processes 
       ! all edge quadrature point evaluations have been completed 
       ! for edges in this subdomain and can be passed to the neighbors 
       ! Send will overlap with internal edge numerical flux calculations
       
       CALL message_send()
       
#else

        WRITE(195,"(ES24.17)") tstage
!         WRITE(195,"(6(ES24.17))") (Qxi(check_iedge,pt)%ptr, pt = 1,nqpte(1))
!         WRITE(195,"(6(ES24.17))") (Qxe(check_iedge,pt)%ptr, pt = 1,nqpte(1))
        WRITE(195,"(6(ES24.17))") (inx(check_iedge,pt), pt = 1,nqpte(1)) 
   
#endif     
   
                  
        
ed_points2: DO pt = 1,nqpte(1) ! Compute numerical fluxes for all edges
              
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)
                
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr*icfac(ed,pt)), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr*icfac(ed,pt)))
              ENDDO
                                        
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)            
!                 Hhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))        
                Zhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(Ze(ed,pt)%ptr - Zi(ed,pt)%ptr))           
              ENDDO
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)       
                Qxhatv(ed) = .5d0*(inx(ed,pt)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))                                     
              ENDDO
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)
                Qyhatv(ed) = .5d0*(inx(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed,pt)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))                                     
              ENDDO
!DIR$ IVDEP              
              DO ed = nfblk(1,blk),nfblk(2,blk) 
!                 Hfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  detJe_in(ed,pt)*Hhatv(ed)
                Zfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Zhatv(ed)
                Zfi(ed,pt)%ptr =  detJe_in(ed,pt)*Zhatv(ed)                
              ENDDO          
!DIR$ IVDEP              
              DO ed = nfblk(1,blk),nfblk(2,blk)                              
                Qxfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP              
              DO ed = nfblk(1,blk),nfblk(2,blk)                                    
                Qyfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qyhatv(ed)
              ENDDO               

        ENDDO ed_points2      
      ENDDO              
   
!$OMP end do 
            
!$OMP do         
      DO blk = 1,nrblk
       DO pt = 1,nqpte(1)
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr*icfac(ed,pt)), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr*icfac(ed,pt)))
              ENDDO             
        
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
!                 Hhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
!                                         - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
                Zhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(Ze(ed,pt)%ptr - Zi(ed,pt)%ptr))           
              ENDDO                                        
     
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
                Qxhatv(ed) = .5d0*(inx(ed,pt)*(xmi(ed,pt)%ptr + xme(ed,pt)%ptr) + iny(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr)  &
                                        - const(ed)*(Qxe(ed,pt)%ptr - Qxi(ed,pt)%ptr))                                 
              ENDDO                                        
            
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
                Qyhatv(ed) = .5d0*(inx(ed,pt)*(xymi(ed,pt)%ptr + xyme(ed,pt)%ptr) + iny(ed,pt)*(ymi(ed,pt)%ptr + yme(ed,pt)%ptr)  &
                                        - const(ed)*(Qye(ed,pt)%ptr - Qyi(ed,pt)%ptr))                                       
              ENDDO                

!DIR$ IVDEP              
              DO ed = rnfblk(1,blk),rnfblk(2,blk)                                     
!                 Hfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Hhatv(ed)
!                 Hfi(ed,pt)%ptr =  detJe_in(ed,pt)*Hhatv(ed)   
                Zfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Zhatv(ed)
                Zfi(ed,pt)%ptr =  detJe_in(ed,pt)*Zhatv(ed)                 
              ENDDO   
!DIR$ IVDEP                              
              DO ed = rnfblk(1,blk),rnfblk(2,blk)                                          
                Qxfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qxhatv(ed)
                Qxfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qxhatv(ed)
              ENDDO   
!DIR$ IVDEP              
              DO ed = rnfblk(1,blk),rnfblk(2,blk)                                             
                Qyfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Qyhatv(ed)
                Qyfi(ed,pt)%ptr =  detJe_in(ed,pt)*Qyhatv(ed)
              ENDDO

        ENDDO            
      ENDDO        
        
!$OMP end do        
        
#ifdef CMPI      

      CALL MPI_WAITALL(2*nproc_sr,solreq,MPI_STATUSES_IGNORE,ierr)
      
      WRITE(195,"(ES24.17)") tstage
!       WRITE(195,"(6(ES24.17))") (Qxri(match_iedge,pt)%ptr, pt = 1,nqpte(1))
!       WRITE(195,"(6(ES24.17))") (Qxre(match_iedge,pt)%ptr, pt = 1,nqpte(1))
      WRITE(195,"(6(ES24.17))") (rnx(match_iedge,pt), pt = 1,nqpte(1))
      
      DO pt = 1,nqpte(1)
      
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
!           Hhatv(ed) = .5d0*(rnx(ed,pt)*(Qxri(ed,pt)%ptr + Qxre(ed,pt)%ptr) + rny(ed,pt)*(Qyri(ed,pt)%ptr + Qyre(ed,pt)%ptr) &
!                                         - const(ed)*(Hre(ed,pt)%ptr - Hri(ed,pt)%ptr))        
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
!           Hfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Hhatv(ed)
          Zfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Zhatv(ed)          
        ENDDO   
        DO ed = 1,nred                                         
          Qxfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qxhatv(ed)        
        ENDDO        
        DO ed = 1,nred 
          Qyfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qyhatv(ed)        
        ENDDO
      ENDDO
#endif              




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

!               Hflux(el_in,gp_in) = detJe(ged,pt)*Hhat
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
                ind = (pt-1)*nfbfr + bfr
                arg = fbfreq(bfr)*(tstage - real(INT(tstage/fbper(bfr)),pres)*fbper(bfr)) + fbeq(bfr)
                Qn = Qn + fbamp_qpt(ind,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(ind,ed))
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
                ind = (pt-1)*nobfr + bfr
                arg = obfreq(bfr)*(tstage-real(INT(tstage/obper(bfr)),pres)*obper(bfr)) + obeq(bfr)
                Z_ex = Z_ex + obamp_qpt(ind,ed)*obnfact(bfr)*ramp*COS(arg-obph_qpt(ind,ed))
              ENDDO
!               H_ex = H_ex + obdepth_qpt(ed,pt)

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




!$OMP do
      DO blk = 1,npart   
        DO et = 1,nel_type
        
          IF (et == 3) THEN 
            ete = 1
          ELSE IF (et == 4) THEN
            ete = 2
          ELSE   
            ete = et
          ENDIF
          
          
          IF (npartet(et,blk) > 0) THEN    
          
            DO pt = 1,nverts(et)*nqpte(ete)
              DO l = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
                DO el = elblk(1,blk,et),elblk(2,blk,et)
                  rhsZ(el,l)  = rhsZ(el,l)  - Zflux(el,pt)*phie_int(l,pt,ete)
                  rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt,ete)
                  rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt,ete)                   
                ENDDO
              ENDDO                                    
            ENDDO  
 
            SELECT CASE(et)
               
              CASE(1)
                DO l = 1,ndof(et)
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    MirhsZ(el,l)  = MirhsZ(el,l)  + mmi(el,1)*rhsZ(el,l) 
                    MirhsQx(el,l) = MirhsQx(el,l) + mmi(el,1)*rhsQx(el,l) 
                    MirhsQy(el,l) = MirhsQy(el,l) + mmi(el,1)*rhsQy(el,l) 
                  ENDDO
                ENDDO
 
              CASE DEFAULT
              
                m = 1
                DO i = 1,ndof(et)
                  DO j = 1,ndof(et)
                    DO el = elblk(1,blk,et),elblk(2,blk,et)
                      MirhsZ(el,i)  = MirhsZ(el,i)  + mmi(el,m)*rhsZ(el,j) 
                      MirhsQx(el,i) = MirhsQx(el,i) + mmi(el,m)*rhsQx(el,j) 
                      MirhsQy(el,i) = MirhsQy(el,i) + mmi(el,m)*rhsQy(el,j) 
                    ENDDO
                    m = m + 1
                  ENDDO
                ENDDO            
            END SELECT
            
          ENDIF
        ENDDO
      ENDDO       
!$OMP enddo



!$OMP end parallel         

      RETURN
      END SUBROUTINE rhs2
