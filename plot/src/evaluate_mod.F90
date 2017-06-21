      MODULE evaluate_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type

      IMPLICIT NONE


      CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_solution(el,et,sol_type,snap,sol_val,npts,r,s,phi_sol,plim)

      USE plot_globals, ONLY: hb,Z,Qx,Qy,el_size,g
      USE read_dginp, ONLY: p,hbp,h0
      USE basis, ONLY: element_basis,linear_basis,dgswem_basis        
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: sol_type
      INTEGER, INTENT(IN) :: snap
      REAL(rp), DIMENSION(:), INTENT(OUT) :: sol_val
      INTEGER, INTENT(IN) :: npts
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: r
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: s   
      REAL(rp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: phi_sol   
      INTEGER, INTENT(IN), OPTIONAL :: plim

      INTEGER :: pt,dof
      INTEGER :: mndf,ndfb,ndfs
      REAL(rp) :: H,u_vel,v_vel
      REAL(rp) :: hb_val(npts),zeta_val(npts)
      REAL(rp) :: Qx_val(npts),Qy_val(npts)
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phis
      REAL(rp) :: depth_min

      mndf = (p+1)**2
      

      
      
      
      IF (PRESENT(r) .and. PRESENT(s)) THEN
       
          ALLOCATE(phis(mndf,npts))           
#ifdef adcirc
          CALL linear_basis(ndfs,npts,r,s,phis) 
#elif dgswem
          CALL dgswem_basis(p,ndfs,npts,r,s,phis)  
#else 
          CALL element_basis(et,p,ndfs,npts,r,s,phis)          
#endif          
      
      ENDIF

      
      IF (PRESENT(phi_sol)) THEN
       
        ALLOCATE(phis(mndf,npts))  
        DO pt = 1,npts
          DO dof = 1,mndf
            phis(dof,pt) = phi_sol(dof,pt)
          ENDDO
        ENDDO     
      
      ENDIF
      
      
      IF (PRESENT(plim)) THEN
        IF (mod(et,2) == 1) THEN
          ndfs = (plim+1)*(plim+2)/2 
          ndfb = min(ndfs,(hbp+1)*(hbp+2)/2)
        ELSE IF (mod(et,2) == 0) THEn
          ndfs = (plim+1)**2
          ndfb = min(ndfs,(hbp+1)**2)
        ENDIF            
      ELSE      
        IF (mod(et,2) == 1) THEN
          ndfs = (p+1)*(p+2)/2 
          ndfb = (hbp+1)*(hbp+2)/2
        ELSE IF (mod(et,2) == 0) THEn
          ndfs = (p+1)**2
          ndfb = (hbp+1)**2
        ENDIF
      ENDIF      
      
      
      
      
      
      IF (sol_type == 2 .or. sol_type == 4 .or. sol_type == 5) THEN                
        
        DO pt = 1,npts
          hb_val(pt) = 0d0
          DO dof = 1,ndfb
            hb_val(pt) = hb_val(pt) + phis(dof,pt)*hb(dof,el,1)
          ENDDO
        ENDDO
            
        IF (sol_type == 2) THEN
          depth_min = 1d10
          DO pt = 1,npts 
            sol_val(pt) = hb_val(pt)
            
            IF (sol_val(pt) < depth_min) THEN
              depth_min = sol_val(pt)
            ENDIF
          ENDDO
          
          IF (depth_min < h0) THEN
            PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb min = ", depth_min      
          ENDIF           
        ENDIF
      ENDIF                         
     
      IF (sol_type == 3 .or. sol_type == 4) THEN                
          
        DO pt = 1,npts
          zeta_val(pt) = 0d0
          DO dof = 1,ndfs
            zeta_val(pt) = zeta_val(pt) + phis(dof,pt)*Z(dof,el,snap)
          ENDDO
        ENDDO
            
        IF (sol_type == 3) THEN
          DO pt = 1,npts
            sol_val(pt) = zeta_val(pt)
          ENDDO
        ENDIF
      ENDIF
          
      IF (sol_type == 4) THEN
        DO pt = 1,npts
          Qx_val(pt) = 0d0
          Qy_val(pt) = 0d0
          DO dof = 1,ndfs
            Qx_val(pt) = Qx_val(pt) + phis(dof,pt)*Qx(dof,el,snap)
            Qy_val(pt) = Qy_val(pt) + phis(dof,pt)*Qy(dof,el,snap)                
          ENDDO
        ENDDO
          
        DO pt = 1,npts 
          H = zeta_val(pt)+hb_val(pt)
          u_vel = Qx_val(pt)/H
          v_vel = Qy_val(pt)/H
          sol_val(pt) = sqrt(u_vel**2+v_vel**2)
        ENDDO
      ENDIF                

      IF (sol_type == 5) THEN

        DO pt = 1,npts 
          sol_val(pt) = el_size(el)/sqrt(g*hb_val(pt))                    
          sol_val(pt) = sol_val(pt)/(2d0*real(p,rp)+1d0)
        ENDDO        

      ENDIF
      
          
      
      RETURN
      END SUBROUTINE evaluate_solution   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
!       SUBROUTINE evaluate_solution(el,et,sol_type,snap,sol_val,npts,r,s,phi_hb,phi_sol,ndf_hb,ndf_sol)
! 
!       USE plot_globals, ONLY: hb,Z,Qx,Qy
!       USE read_dginp, ONLY: p,hbp,h0
!       USE basis, ONLY: element_basis,linear_basis,dgswem_basis        
!       
!       IMPLICIT NONE
!       
!       INTEGER, INTENT(IN) :: el
!       INTEGER, INTENT(IN) :: et
!       INTEGER, INTENT(IN) :: sol_type
!       INTEGER, INTENT(IN) :: snap
!       REAL(rp), DIMENSION(:), INTENT(OUT) :: sol_val
!       INTEGER, INTENT(IN) :: npts
!       REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: r
!       REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: s   
!       REAL(rp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: phi_hb
!       REAL(rp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: phi_sol   
!       INTEGER, INTENT(IN), OPTIONAL :: ndf_hb
!       INTEGER, INTENT(IN), OPTIONAL :: ndf_sol
! 
!       INTEGER :: pt,dof
!       INTEGER :: mndf,ndfb,ndfs
!       INTEGER :: dofmap(6),ndif_dof      
!       REAL(rp) :: H,u_vel,v_vel
!       REAL(rp) :: hb_val(npts),zeta_val(npts)
!       REAL(rp) :: Qx_val(npts),Qy_val(npts)
!       REAL(rp) :: hb_lo(npts),zeta_lo(npts)
!       REAL(rp) :: Qx_lo(npts),Qy_lo(npts)      
!       REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phib
!       REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phis
!       REAL(rp) :: depth_min
!       
! !       ndif_dof = 3
! !       dofmap = (/ 1, 2, 5, 0, 0, 0 /)
!       ndif_dof = 6
!       dofmap = (/ 1, 2, 3, 5, 6, 8 /)
! 
!       mndf = (p+1)**2
!       
!       IF (PRESENT(r) .and. PRESENT(s)) THEN
!       
!         IF (sol_type == 2 .or. sol_type == 4) THEN   
!           ALLOCATE(phib(mndf,npts))  
! #ifdef adcirc
!           CALL linear_basis(ndfb,npts,r,s,phib) 
! #elif dgswem
!           CALL dgswem_basis(hbp,ndfb,npts,r,s,phib) 
! #else
!           CALL element_basis(et,hbp,ndfb,npts,r,s,phib)        
! #endif          
!         ENDIF
!       
!         IF (sol_type == 3 .or. sol_type == 4) THEN        
!           ALLOCATE(phis(mndf,npts))           
! #ifdef adcirc
!           CALL linear_basis(ndfs,npts,r,s,phis) 
! #elif dgswem
!           CALL dgswem_basis(p,ndfs,npts,r,s,phis)  
! #else 
!           CALL element_basis(et,p,ndfs,npts,r,s,phis)          
! #endif          
!         ENDIF
!       
!       ENDIF
!       
!       
!       IF (PRESENT(phi_hb)) THEN
!         
!         IF (sol_type == 2 .or. sol_type == 4) THEN   
!           ALLOCATE(phib(mndf,npts))    
!           DO pt = 1,npts
!             DO dof = 1,mndf
!               phib(dof,pt) = phi_hb(dof,pt)
!             ENDDO
!           ENDDO
!           
!           ndfb = ndf_hb
!         ENDIF
!          
!       ENDIF
!       
!       IF (PRESENT(phi_sol)) THEN
!       
!         IF (sol_type == 3 .or. sol_type == 4) THEN        
!           ALLOCATE(phis(mndf,npts))  
!           DO pt = 1,npts
!             DO dof = 1,mndf
!               phis(dof,pt) = phi_sol(dof,pt)
!             ENDDO
!           ENDDO  
!           
!           ndfs = ndf_sol
!         ENDIF      
!       
!       ENDIF
!       
!       
!       
!       
!       
!       IF (sol_type == 2 .or. sol_type == 4) THEN                
!         
!         DO pt = 1,npts
!           hb_val(pt) = 0d0
!           DO dof = 1,ndfb
!             hb_val(pt) = hb_val(pt) + phib(dof,pt)*hb(dof,el,1)
!           ENDDO
!           
!           !!!!!!!!!!!!!!!!!!!!!!!!!!
!           hb_lo(pt) = 0d0
!           DO dof = 1,ndif_dof
!             hb_lo(pt) = hb_lo(pt) + phib(dofmap(dof),pt)*hb(dofmap(dof),el,1)
!           ENDDO
! !           DO dof = 1,ndif_dof
! !             hb_val(pt) = hb_val(pt) - phib(dofmap(dof),pt)*hb(dofmap(dof),el,1)
! !           ENDDO
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ENDDO
!             
!         IF (sol_type == 2) THEN
!           depth_min = 1d10
!           DO pt = 1,npts 
! !             sol_val(pt) = hb_val(pt)
!             sol_val(pt) = abs(hb_val(pt) - hb_lo(pt))
! !             sol_val(pt) = hb_lo(pt)            
!             
!             IF (sol_val(pt) < depth_min) THEN
!               depth_min = sol_val(pt)
!             ENDIF
!           ENDDO
!           
!           IF (depth_min < h0) THEN
!             PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb min = ", depth_min      
!           ENDIF           
!         ENDIF
!       ENDIF                         
!      
!       IF (sol_type == 3 .or. sol_type == 4) THEN                
!           
!         DO pt = 1,npts
!           zeta_val(pt) = 0d0
!           DO dof = 1,ndfs
!             zeta_val(pt) = zeta_val(pt) + phis(dof,pt)*Z(dof,el,snap)
!           ENDDO
!           
!           !!!!!!!!!!!!!!!!!!!!!!!!!
!           zeta_lo(pt) = 0d0
!           DO dof = 1,ndif_dof
!             zeta_lo(pt) = zeta_lo(pt) + phis(dofmap(dof),pt)*Z(dofmap(dof),el,snap)
!           ENDDO
! !           DO dof = 1,ndif_dof
! !             zeta_val(pt) = zeta_val(pt) - phis(dofmap(dof),pt)*Z(dofmap(dof),el,snap)
! !           ENDDO
!           !!!!!!!!!!!!!!!!!!!!!!!!!
!         ENDDO
!             
!         IF (sol_type == 3) THEN
!           DO pt = 1,npts
! !             sol_val(pt) = zeta_val(pt)
!             sol_val(pt) = zeta_val(pt) - zeta_lo(pt)
! !             sol_val(pt) = zeta_lo(pt)            
!           ENDDO
!         ENDIF
!       ENDIF
!           
!       IF (sol_type == 4) THEN
!         DO pt = 1,npts
!           Qx_val(pt) = 0d0
!           Qy_val(pt) = 0d0
!           DO dof = 1,ndfs
!             Qx_val(pt) = Qx_val(pt) + phis(dof,pt)*Qx(dof,el,snap)
!             Qy_val(pt) = Qy_val(pt) + phis(dof,pt)*Qy(dof,el,snap)                
!           ENDDO
!           
!           !!!!!!!!!!!!!!!!!!!!!
!           Qx_lo(pt) = 0d0
!           Qy_lo(pt) = 0d0
!           DO dof = 1,ndif_dof
!             Qx_lo(pt) = Qx_lo(pt) + phis(dofmap(dof),pt)*Qx(dofmap(dof),el,snap)
!             Qy_lo(pt) = Qy_lo(pt) + phis(dofmap(dof),pt)*Qy(dofmap(dof),el,snap)
!           ENDDO
! 
! !           DO dof = 1,ndif_dof
! !             Qx_val(pt) = Qx_val(pt) - phis(dofmap(dof),pt)*Qx(dofmap(dof),el,snap)
! !             Qy_val(pt) = Qy_val(pt) - phis(dofmap(dof),pt)*Qy(dofmap(dof),el,snap)
! !           ENDDO
!           !!!!!!!!!!!!!!!!!!!!!
!         ENDDO
!           
!         DO pt = 1,npts 
!           H = zeta_val(pt)+hb_val(pt)
!           u_vel = Qx_val(pt)/H
!           v_vel = Qy_val(pt)/H
!           sol_val(pt) = sqrt(u_vel**2+v_vel**2)
!           
!           H = zeta_lo(pt)+hb_lo(pt)
!           u_vel = Qx_lo(pt)/H
!           v_vel = Qy_lo(pt)/H
!           sol_val(pt) = abs(sol_val(pt) - sqrt(u_vel**2+v_vel**2))
! !           sol_val(pt) = abs(sol_val(pt) - sqrt(u_vel**2+v_vel**2))/sqrt(u_vel**2+v_vel**2)           
! !           sol_val(pt) = (sol_val(pt) - sqrt(u_vel**2+v_vel**2))/sol_val(pt)              
!         ENDDO
!       ENDIF                
! 
!       
!       RETURN
!       END SUBROUTINE evaluate_solution        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE evaluate_basis(p,nord,mnpp,mndof,nel_type,npplt,r,s,ndof,phi)
      
      USE basis, ONLY: element_basis,linear_basis,dgswem_basis
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(IN) :: nord
      INTEGER, INTENT(IN) :: mnpp
      INTEGER, INTENT(IN) :: mndof
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: s
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ndof
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: phi
      
      INTEGER :: ord
      INTEGER :: i,j
      INTEGER :: et
      INTEGER :: n

      
      ALLOCATE(phi(mndof,mnpp,nel_type*nord))
      ALLOCATE(ndof(nel_type))

      DO ord = 1,nord
        DO et = 1,nel_type
          i = (et-1)*nord+ord
          
#ifdef adcirc            
          CALL linear_basis(ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i))             
#elif dgswem
          CALL dgswem_basis(p,ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i))          
#else        
          CALL element_basis(et,p,ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i)) 
#endif
        
        ENDDO          
      ENDDO
      
      RETURN
      END SUBROUTINE evaluate_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      END MODULE evaluate_mod