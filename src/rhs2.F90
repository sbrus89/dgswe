      SUBROUTINE rhs2()

      USE globals,  ONLY: pres,ndof,nel_type, &
                          elblk,nfblk,rnfblk,nrblk,npartet,nverts, &                        
                          Hi,He,Qxi,Qxe,Qyi,Qye,hbqpted, &
                          nqpta,nqpte, &                     
                          check_iedge,check_gedge, &                     
                          const,inx,iny,icfac,detJe_in,detJe_ex,detJe, &
                          nx_pt,ny_pt,Spe, &
                          rhsZ,rhsQx,rhsQy, &
                          MirhsZ,MirhsQx,MirhsQy

                          
      USE read_dginp, ONLY: npart                          
                          
#ifdef CMPI                          
      USE messenger2, ONLY: myrank,message_recieve,message_send, &
                            nred,nproc_sr,match_edge, &
                            solreq,solreq_send,solreq_recv,ierr, &
                            Zri,Zre,Qxri,Qxre,Qyri,Qyre, &
                            rnx,rny,rcfac,detJe_recv,hbr         
                            
      USE mpi                            
#endif       
!$    USE omp_lib                          
                          
                   
      IMPLICIT NONE
      INTEGER :: et,blk
      INTEGER :: ete

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     c Area Integrals
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

!!DIR$ VECTOR ALIGNED
      rhsZ(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQx(:,:) = 0d0
!!DIR$ VECTOR ALIGNED
      rhsQy(:,:) = 0d0   
      
!!DIR$ VECTOR ALIGNED
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
            CALL area_integration(et,elblk(1,blk,et),elblk(2,blk,et),ndof(et),nqpta(et))                        
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
            CALL interior_edge_eval(ete,elblk(1,blk,et),elblk(2,blk,et),ndof(et),nverts(et)*nqpte(ete))     
         ENDIF
        ENDDO
        
#ifdef CMPI    
       ! Post an non-blocking send to all processes 
       ! all edge quadrature point evaluations have been completed 
       ! for edges in this subdomain and can be passed to the neighbors 
       ! Send will overlap with internal edge numerical flux calculations
       
       CALL message_send()
       
#else

!         WRITE(195,"(ES24.17)") tstage
! !         WRITE(195,"(6(ES24.17,1x))") (Qxi(check_iedge,pt)%ptr, pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (Qxe(check_iedge,pt)%ptr, pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (Qyi(check_iedge,pt)%ptr, pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (Qye(check_iedge,pt)%ptr, pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (Zi(check_iedge,pt)%ptr, pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (Ze(check_iedge,pt)%ptr, pt = 1,nqpte(1))
!           WRITE(195,"(6(ES24.17,1x))") (inx(check_iedge,pt), pt = 1,nqpte(1)) 
! !         WRITE(195,"(6(ES24.17,1x))") (iny(check_iedge,pt), pt = 1,nqpte(1)) 
! !         WRITE(195,"(6(ES24.17,1x))") (detJe_ex(check_iedge,pt), pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (detJe_in(check_iedge,pt), pt = 1,nqpte(1))
! !         WRITE(195,"(6(ES24.17,1x))") (hbqpted(check_gedge,pt), pt = 1,nqpte(1))

   
#endif     
   
                  
        CALL interior_edge_nflux(nfblk(1,blk),nfblk(2,blk),nqpte(1))

      ENDDO              
   
!$OMP end do 
            
!$OMP do         
      DO blk = 1,nrblk
        CALL interior_edge_nflux(rnfblk(1,blk),rnfblk(2,blk),nqpte(1))         
      ENDDO        
        
!$OMP end do        
        
#ifdef CMPI      

      CALL MPI_WAITALL(2*nproc_sr,solreq,MPI_STATUSES_IGNORE,ierr)
      
!       WRITE(195,"(ES24.17)") tstage
! !       WRITE(195,"(6(ES24.17,1x))") (Qxri(match_edge,pt)%ptr, pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (Qxre(match_edge,pt)%ptr, pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (Qyri(match_edge,pt)%ptr, pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (Qyre(match_edge,pt)%ptr, pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (Zri(match_edge,pt)%ptr, pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (Zre(match_edge,pt)%ptr, pt = 1,nqpte(1))
!         WRITE(195,"(6(ES24.17,1x))") (rnx(match_edge,pt), pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (rny(match_edge,pt), pt = 1,nqpte(1))
! !      WRITE(195,"(6(ES24.17,1x))") (detJe_recv(match_edge,pt), pt = 1,nqpte(1))
! !       WRITE(195,"(6(ES24.17,1x))") (hbr(match_edge,pt), pt = 1,nqpte(1))

      CALL recieve_edge_nflux(nred,nqpte(1))      

#endif              




      CALL boundary_edge_land()
        
      CALL boundary_edge_flow()        

      CALL boundary_edge_elev()
      




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
          
            CALL edge_integration(ete,elblk(1,blk,et),elblk(2,blk,et),ndof(et),nverts(et)*nqpte(ete))  
 
            CALL linear_solve(et,elblk(1,blk,et),elblk(2,blk,et),ndof(et))
            
          ENDIF
        ENDDO
      ENDDO       
!$OMP enddo



!$OMP end parallel         

      RETURN
      END SUBROUTINE rhs2
