      SUBROUTINE rhs2()

      USE globals,  ONLY: pres,l,ind,el,ne,pt,led,ed,dof,ndof,nel_type, &
                          blk,elblk,nfblk,rnfblk,nrblk,npart,npartet,nverts, &
                          g,pt5g,tstage,ramp,arg,bfr, & 
                          H,Qx,Qy, &
                          rhsH,rhsQx,rhsQy, &
                          Hqpt,Qxqpt,Qyqpt, &
                          nqpta,wpta,phia,phia_int,dpdx,dpdy,nqpte,wpte,phie,phie_int,mmi, &
                          press,recipH,xmom,ymom,xymom,cf,tau,u,v,src_x,src_y,dhbdx,dhbdy, &
                          el_in,el_ex,led_in,led_ex,gp_in,gp_ex, &
                          nx,ny,nx2,ny2,nxny,tx,ty, &
                          H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex,Qn,Qt, &
                          xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
                          Hhat,Qxhat,Qyhat, &
                          Hflux,Qxflux,Qyflux, &
                          nied,iedn,nobed,obedn,nnfbed,nfbedn,nfbed,fbedn, &
                          nobfr,obfreq,obper,obnfact,obamp_qpt,obph_qpt,obeq,obdepth_qpt, &
                          nfbfr,fbfreq,fbper,fbnfact,fbamp_qpt,fbph_qpt,fbeq, &
                          ged2led,ged2el,gel2ael,ged, &
                          pressa,recipHa, & 
                          Hi,He,Qxi,Qxe,Qyi,Qye, &
                          xmi,xme,ymi,yme,xymi,xyme, &
                          Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                          const,inx,iny,detJe_in,detJe_ex,detJe,nx_pt,ny_pt, &
                          Hhatv,Qxhatv,Qyhatv, &
                          MirhsH,MirhsQx,MirhsQy
                          
#ifdef CMPI                          
      USE messenger2, ONLY: message_recieve,message_send, &
                            nred,nproc_sr, &
                            solreq,solreq_send,solreq_recv,ierr, &
                            Hri,Hre,Qxri,Qxre,Qyri,Qyre, &
                            xmri,ymri,xymri,xmre,ymre,xymre, &
                            Hfri,Qxfri,Qyfri, &
                            rnx,rny,detJe_recv         
                            
      USE mpi                            
#endif       
!$    USE omp_lib                          
                          
                   
      IMPLICIT NONE
      INTEGER :: i,j,et,m

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

!!DIR$ VECTOR ALIGNED
      rhsH(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQy(:,:) = 0d0   
      
!!DIR$ VECTOR ALIGNED
      MirhsH(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      MirhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      MirhsQy(:,:) = 0d0         

      
!$OMP parallel default(none)  &
!$OMP             private(ind,blk,pt,dof,el,l,ed,ged,led_in,gp_in,el_in, &
!$OMP                     nx,ny,tx,ty,nx2,ny2,nxny,H_in,H_ex,Qx_in,Qx_ex,Qy_in,Qy_ex, &
!$OMP                     press,recipH,xmom_in,xmom_ex,ymom_in,ymom_ex,xymom_in,xymom_ex, &
!$OMP                     Hhat,Qxhat,Qyhat, &
!$OMP                     Qn,Qt,arg,bfr) &
!$OMP             shared(nqpta,nqpte,ndof,H,Qx,Qy,Hqpt,Qxqpt,Qyqpt, & 
!$OMP                    phia,phie,recipHa,xmom,ymom,xymom,tau,src_x,src_y, &
!$OMP                    dhbdx,dhbdy,cf,phia_int,phie_int,dpdx,dpdy,rhsH,rhsQx,rhsQy, &
!$OMP                    const,Qxi,Qxe,Qyi,Qye,Hi,He,xmi,xme,ymi,yme,xymi,xyme, &
!$OMP                    Hhatv,Qxhatv,Qyhatv,Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
!$OMP                    Hflux,Qxflux,Qyflux,inx,iny,len_area_in,len_area_ex, &
!$OMP                    nnfbed,nfbedn,nfbed,fbedn,nobed,obedn, &
!$OMP                    nfbfr,fbfreq,fbper,fbeq,fbamp_qpt,fbnfact,fbph_qpt, &
!$OMP                    nobfr,obfreq,obper,obeq,obamp_qpt,obnfact,obph_qpt,obdepth_qpt, &
!$OMP                    ramp,tstage, &
!$OMP                    edlen_area,normal, &
!$OMP                    ged2led,gel2ael,ged2el, &
!$OMP                    nblk,npart,nrblk,elblk,edblk,nfblk,rnfblk)


!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

!$OMP do
   

      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
            CALL area_integration(et,elblk(1,blk,et),elblk(2,blk,et),ndof(et),nqpta(et))
          ENDIF
        ENDDO
      ENDDO       

      
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
          IF (npartet(et,blk) > 0) THEN       
ed_points: DO pt = 1,nverts(et)*nqpte(et)

!!DIR$ VECTOR ALIGNED               
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Hqpt(el,pt) = H(el,1)
                Qxqpt(el,pt) = Qx(el,1)
                Qyqpt(el,pt) = Qy(el,1)
              ENDDO

    ed_basis: DO dof = 2,ndof(et)
!!DIR$ VECTOR ALIGNED
                DO el = elblk(1,blk,et),elblk(2,blk,et)   ! Compute solutions at edge quadrature points                
                  Hqpt(el,pt) = Hqpt(el,pt) + H(el,dof)*phie(dof,pt,et)            
                  Qxqpt(el,pt) = Qxqpt(el,pt) + Qx(el,dof)*phie(dof,pt,et)            
                  Qyqpt(el,pt) = Qyqpt(el,pt) + Qy(el,dof)*phie(dof,pt,et)            
                ENDDO

              ENDDO ed_basis

!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)  ! Compute momentum terms       
                recipHa(el) = 1d0/Hqpt(el,pt)
                
                xmom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qxqpt(el,pt)*Qxqpt(el,pt)*recipHa(el) 
                ymom(el,pt) = pt5g*Hqpt(el,pt)*Hqpt(el,pt) + Qyqpt(el,pt)*Qyqpt(el,pt)*recipHa(el) 
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
#endif           
                  
        
ed_points2: DO pt = 1,nqpte(1) ! Compute numerical fluxes for all edges
              
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)
                
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO
                                        
              
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = nfblk(1,blk),nfblk(2,blk)            
                Hhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
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
                Hfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  detJe_in(ed,pt)*Hhatv(ed)
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
                const(ed) = max(abs(Qxi(ed,pt)%ptr*inx(ed,pt) + Qyi(ed,pt)%ptr*iny(ed,pt))/Hi(ed,pt)%ptr + sqrt(g*Hi(ed,pt)%ptr), &
                                abs(Qxe(ed,pt)%ptr*inx(ed,pt) + Qye(ed,pt)%ptr*iny(ed,pt))/He(ed,pt)%ptr + sqrt(g*He(ed,pt)%ptr))
              ENDDO             
        
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO ed = rnfblk(1,blk),rnfblk(2,blk)
                Hhatv(ed) = .5d0*(inx(ed,pt)*(Qxi(ed,pt)%ptr + Qxe(ed,pt)%ptr) + iny(ed,pt)*(Qyi(ed,pt)%ptr + Qye(ed,pt)%ptr) &
                                        - const(ed)*(He(ed,pt)%ptr - Hi(ed,pt)%ptr))
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
                Hfe(ed,pt)%ptr = -detJe_ex(ed,pt)*Hhatv(ed)
                Hfi(ed,pt)%ptr =  detJe_in(ed,pt)*Hhatv(ed)               
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
        
#ifdef CMPI      

      CALL MPI_WAITALL(2*nproc_sr,solreq,MPI_STATUSES_IGNORE,ierr)

      
      DO pt = 1,nqpte(1)
      
!!DIR$ VECTOR ALIGNED      
        DO ed = 1,nred
          const(ed) = max(abs(Qxri(ed,pt)%ptr*rnx(ed,pt) + Qyri(ed,pt)%ptr*rny(ed,pt))/Hri(ed,pt)%ptr + sqrt(g*Hri(ed,pt)%ptr), &
                          abs(Qxre(ed,pt)%ptr*rnx(ed,pt) + Qyre(ed,pt)%ptr*rny(ed,pt))/Hre(ed,pt)%ptr + sqrt(g*Hre(ed,pt)%ptr))          
        ENDDO
        

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
        DO ed = 1,nred
          Hhatv(ed) = .5d0*(rnx(ed,pt)*(Qxri(ed,pt)%ptr + Qxre(ed,pt)%ptr) + rny(ed,pt)*(Qyri(ed,pt)%ptr + Qyre(ed,pt)%ptr) &
                                        - const(ed)*(Hre(ed,pt)%ptr - Hri(ed,pt)%ptr))
        ENDDO

!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED        
        DO ed = 1,nred
          recipHa(ed) = 1d0/Hre(ed,pt)%ptr
          
          xmre(ed) = pt5g*Hre(ed,pt)%ptr*Hre(ed,pt)%ptr + Qxre(ed,pt)%ptr*Qxre(ed,pt)%ptr*recipHa(ed)
          ymre(ed) = pt5g*Hre(ed,pt)%ptr*Hre(ed,pt)%ptr + Qyre(ed,pt)%ptr*Qyre(ed,pt)%ptr*recipHa(ed)
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
          Hfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Hhatv(ed)
        ENDDO   
        DO ed = 1,nred                                         
          Qxfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qxhatv(ed)        
        ENDDO        
        DO ed = 1,nred 
          Qyfri(ed,pt)%ptr =  detJe_recv(ed,pt)*Qyhatv(ed)        
        ENDDO
      ENDDO
#endif              
        
      ENDDO
!$OMP end do


!$OMP do     
            ! No normal flow boundary condition 
            DO ed = 1,nnfbed
              DO pt = 1,nqpte(3)

              ged = nfbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte(3) + pt

              el_in = gel2ael(ged2el(1,ged))

              nx = nx_pt(ged,pt)
              ny = ny_pt(ged,pt)

              nx2 = nx*nx
              ny2 = ny*ny
              nxny = nx*ny

              H_in = Hqpt(el_in,gp_in)
              H_ex = H_in

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qx_ex = Qx_in*(ny2-nx2) - 2d0*nxny*Qy_in
              Qy_ex = Qy_in*(nx2-ny2) - 2d0*nxny*Qx_in

              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = detJe(ged,pt)*Hhat

              Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat

              Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
              ENDDO
            ENDDO
!$OMP end do

!$OMP do

            ! Flow specified boundary edges
            DO ed = 1,nfbed
              DO pt = 1,nqpte(1)

              ged = fbedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte(1) + pt

              el_in = gel2ael(ged2el(1,ged))
              
              nx = nx_pt(ged,pt)
              ny = ny_pt(ged,pt)

              tx = -ny
              ty = nx

              H_in = Hqpt(el_in,gp_in)
              H_ex = H_in

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qn = 0d0
              Qt= 0d0
              DO bfr = 1,nfbfr
                ind = (pt-1)*nfbfr + bfr
                arg = fbfreq(bfr)*(tstage - INT(tstage/fbper(bfr))*fbper(bfr)) + fbeq(bfr)
                Qn = Qn + fbamp_qpt(ind,ed)*fbnfact(bfr)*ramp*COS(arg-fbph_qpt(ind,ed))
              ENDDO

              Qn = -Qn

              Qx_ex = ( ty*Qn - ny*Qt)/(nx*ty-ny*tx)
              Qy_ex = (-tx*Qn + nx*Qt)/(nx*ty-ny*tx)
 
              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = detJe(ged,pt)*Hhat
              
              Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat

              Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat
              
              ENDDO
            ENDDO
            
!$OMP end do

!$OMP do

             ! Open boundary edges (elevation specified)
            DO ed = 1,nobed
              DO pt = 1,nqpte(1)

              ged = obedn(ed)

              led_in = ged2led(1,ged)

              gp_in = (led_in-1)*nqpte(1) + pt

              el_in = gel2ael(ged2el(1,ged))

              nx = nx_pt(ged,pt)
              ny = ny_pt(ged,pt)

              H_in = Hqpt(el_in,gp_in)

              H_ex = 0d0
              DO bfr = 1,nobfr
                ind = (pt-1)*nobfr + bfr
                arg = obfreq(bfr)*(tstage-INT(tstage/obper(bfr))*obper(bfr)) + obeq(bfr)
                H_ex = H_ex + obamp_qpt(ind,ed)*obnfact(bfr)*ramp*COS(arg-obph_qpt(ind,ed))
              ENDDO
              H_ex = H_ex + obdepth_qpt(ed,pt)

              Qx_in = Qxqpt(el_in,gp_in)
              Qy_in = Qyqpt(el_in,gp_in)

              Qx_ex = Qx_in
              Qy_ex = Qy_in
 
              CALL numerical_flux(Qx_in,Qy_in,H_in,Qx_ex,Qy_ex,H_ex,nx,ny,Qxhat,Qyhat,Hhat)

              Hflux(el_in,gp_in) = detJe(ged,pt)*Hhat

              Qxflux(el_in,gp_in) = detJe(ged,pt)*Qxhat

              Qyflux(el_in,gp_in) = detJe(ged,pt)*Qyhat                            
              ENDDO
             ENDDO
!$OMP end do        


      DO blk = 1,npart   
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN    
          
            DO pt = 1,nverts(et)*nqpte(et)
              DO l = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
                DO el = elblk(1,blk,et),elblk(2,blk,et)
                  rhsH(el,l)  = rhsH(el,l)  - Hflux(el,pt)*phie_int(l,pt,et)
                  rhsQx(el,l) = rhsQx(el,l) - Qxflux(el,pt)*phie_int(l,pt,et)
                  rhsQy(el,l) = rhsQy(el,l) - Qyflux(el,pt)*phie_int(l,pt,et)                   
                ENDDO
              ENDDO                                    
            ENDDO  
 
            SELECT CASE(et)
               
              CASE(1)
                DO l = 1,ndof(et)
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    MirhsH(el,l)  = MirhsH(el,l)  + mmi(el,1)*rhsH(el,l) 
                    MirhsQx(el,l) = MirhsQx(el,l) + mmi(el,1)*rhsQx(el,l) 
                    MirhsQy(el,l) = MirhsQy(el,l) + mmi(el,1)*rhsQy(el,l) 
                  ENDDO
                ENDDO
 
              CASE DEFAULT
                m = 1
                DO i = 1,ndof(et)
                  DO j = 1,ndof(et)
                    DO el = elblk(1,blk,et),elblk(2,blk,et)
                      MirhsH(el,i)  = MirhsH(el,i)  + mmi(el,m)*rhsH(el,j) 
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




!$OMP end parallel         

      RETURN
      END SUBROUTINE rhs2