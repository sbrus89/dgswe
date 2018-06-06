      MODULE evaluate_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type,solution_type

      IMPLICIT NONE


      CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_solution(el,sol_type,snap,sol,sol_val,npts,r,s,phi_sol,plim)

      USE basis, ONLY: element_basis,linear_basis,linear_quad_basis,dgswem_basis        
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: sol_type
      INTEGER, INTENT(IN) :: snap
      TYPE(solution_type), INTENT(IN) :: sol
      REAL(rp), DIMENSION(:), INTENT(OUT) :: sol_val
      INTEGER, INTENT(IN) :: npts
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: r
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: s   
      REAL(rp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: phi_sol   
      INTEGER, INTENT(IN), OPTIONAL :: plim

      INTEGER :: pt,dof,et
      INTEGER :: mndf,ndfb,ndfs
      INTEGER :: sol_snap
      REAL(rp) :: H,u_vel,v_vel
      REAL(rp) :: hb_val(npts),zeta_val(npts)
      REAL(rp) :: Qx_val(npts),Qy_val(npts)
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phis
      REAL(rp) :: depth_min

      mndf = (sol%p+1)**2
      

      et = sol%el_type(el)
      
      
      IF (PRESENT(r) .and. PRESENT(s)) THEN
       
        ALLOCATE(phis(mndf,npts))           
        IF (sol%output_type == "adcirc") THEN
          IF (mod(et,2) == 1) THEN          
            CALL linear_basis(ndfs,npts,r,s,phis) 
          ELSE IF (mod(et,2) == 0) THEN
            CALL linear_quad_basis(ndfs,npts,r,s,phis)           
          ENDIF
        ELSE IF (sol%output_type == "dgswem") THEN
          CALL dgswem_basis(sol%p,ndfs,npts,r,s,phis)  
        ELSE IF (sol%output_type == "dgswe") THEN
          CALL element_basis(et,sol%p,ndfs,npts,r,s,phis)          
        ENDIF
              
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
          ndfb = min(ndfs,(sol%hbp+1)*(sol%hbp+2)/2)
        ELSE IF (mod(et,2) == 0) THEn
          ndfs = (plim+1)**2
          ndfb = min(ndfs,(sol%hbp+1)**2)
        ENDIF            
      ELSE      
        IF (mod(et,2) == 1) THEN
          ndfs = (sol%p+1)*(sol%p+2)/2 
          ndfb = (sol%hbp+1)*(sol%hbp+2)/2
        ELSE IF (mod(et,2) == 0) THEn
          ndfs = (sol%p+1)**2
          ndfb = (sol%hbp+1)**2
        ENDIF
      ENDIF      
      
      
      ! Account for initial condition in dgswe output      
      IF (sol%output_type == "dgswe") THEN
        sol_snap = snap + 1
      ELSE 
        sol_snap = snap
      ENDIF          
      
      
      
      IF (sol_type == 2 .or. sol_type == 4 .or. sol_type == 5) THEN                
        
        DO pt = 1,npts
          hb_val(pt) = 0d0
          DO dof = 1,ndfb
            hb_val(pt) = hb_val(pt) + phis(dof,pt)*sol%hb(dof,el,1)
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
          
          IF (depth_min < sol%h0) THEN
            PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb min = ", depth_min      
          ENDIF           
        ENDIF
      ENDIF                         
     
      IF (sol_type == 3 .or. sol_type == 4) THEN                
          
        DO pt = 1,npts
          zeta_val(pt) = 0d0
          DO dof = 1,ndfs
            zeta_val(pt) = zeta_val(pt) + phis(dof,pt)*sol%Z(dof,el,sol_snap)
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
            Qx_val(pt) = Qx_val(pt) + phis(dof,pt)*sol%Qx(dof,el,sol_snap)
            Qy_val(pt) = Qy_val(pt) + phis(dof,pt)*sol%Qy(dof,el,sol_snap)                
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
          sol_val(pt) = sol%el_size(el)/sqrt(sol%g*hb_val(pt))                    
          sol_val(pt) = sol_val(pt)/(2d0*real(sol%p,rp)+1d0)
        ENDDO        

      ENDIF
      
          
      
      RETURN
      END SUBROUTINE evaluate_solution   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        


      SUBROUTINE evaluate_basis(output_type,p,nord,mnpp,mndof,nel_type,npplt,r,s,ndof,phi)
      
      USE basis, ONLY: element_basis,linear_basis,linear_quad_basis,dgswem_basis
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: output_type
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
           
          IF (output_type == "adcirc") THEN
            IF (mod(et,2) == 1) THEN
              CALL linear_basis(ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i))             
            ELSE IF (mod(et,2) == 0) THEN
              CALL linear_quad_basis(ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i))            
            ENDIF
          ELSE IF (output_type == "dgswem") THEN
            CALL dgswem_basis(p,ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i))          
          ELSE IF (output_type == "dgswe") THEN
            CALL element_basis(et,p,ndof(et),npplt(i),r(:,i),s(:,i),phi(:,:,i)) 
          ENDIF
        
        ENDDO          
      ENDDO
      
      RETURN
      END SUBROUTINE evaluate_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate_plotting_nodes(nel_type,pplt_low,pplt_skip,nord,mnpp,np,mnnds,pplt,npplt,r,s,psic,nptri,rect)

      USE plot_globals, ONLY: figure_width
      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE triangulation, ONLY: reference_element_delaunay

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, INTENT(IN) :: pplt_low
      INTEGER, INTENT(IN) :: pplt_skip
      INTEGER, INTENT(IN) :: nord
      INTEGER, INTENT(IN) :: mnpp
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      INTEGER, INTENT(IN) :: mnnds
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: pplt
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: npplt
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: r,s
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: psic
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nptri
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: rect

      INTEGER :: ord,et
      INTEGER :: el,nd
      INTEGER :: space
      INTEGER :: i,nnd
      INTEGER :: ncall
      INTEGER :: plot 
      CHARACTER(3) :: nout
      CHARACTER(25) :: fname

      ! Order of plotting nodal set, and corresponding number of nodes
      ALLOCATE(pplt(nel_type*nord),npplt(nel_type*nord)) 
      
      ! Reference element coordinates for plotting nodes
      ALLOCATE(r(mnpp,nel_type*nord),s(mnpp,nel_type*nord))         

      ! Shape functions evaluated at plotting nodes
      ALLOCATE(psic(mnnds,mnpp,nel_type*nord))

      ! Sub-element connectivity table for reference element plotting nodes
      ALLOCATE(rect(3,3*mnpp,nel_type*nord))      

      ! Number of plotting sub-elements 
      ALLOCATE(nptri(nel_type*nord))

      plot = 0
      ncall = 0 
      DO et = 1,nel_type 
              
        DO ord = 1,nord

          i = (et-1)*nord+ord
          pplt(i) = (ord-1)*pplt_skip+pplt_low

          CALL element_nodes(et,space,pplt(i),npplt(i),r(:,i),s(:,i))                  
          CALL shape_functions_area_eval(et,np(et),nnd,npplt(i),r(:,i),s(:,i),psic(:,:,i))  
          CALL reference_element_delaunay(et,pplt(i),npplt(i),r(:,i),s(:,i),nptri(i),rect(:,:,i))        
          
          PRINT("(4(A,I4))"), "  number of additional nodes/sub-triangles: ", npplt(i),"/",nptri(i) 
          

!           ! Plot reference element triangularization
!           IF (plot == 1) THEN
!             ncall = ncall + 1
!             WRITE(nout,"(I3.3)") ncall 
!             fname = "ref_el_"//nout//".ps"
!             CALL plot_ref_el(fname,figure_width,et,np(et),nptri(i),rect(:,:,i),r(:,i),s(:,i))
! 
!             DO el = 1,nptri(i)
!               PRINT "(4(I5))", el,(rect(nd,el,i), nd=1,3)
!             ENDDO
!             PRINT*, "" 
!           ENDIF
          

        ENDDO         
        
      ENDDO                                    
      RETURN
      END SUBROUTINE evaluate_plotting_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE contour_line_newton(sol,type_flag,sol_lev,el,et,sol_snap,r1,r2,s1,s2,re,se)
      
      USE basis, ONLY: element_basis,linear_basis,dgswem_basis      
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(IN) :: sol    
      INTEGER, INTENT(IN) :: type_flag
      REAL(rp), INTENT(IN) :: sol_lev
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: sol_snap
      REAL(rp), INTENT(IN) :: r1,r2
      REAL(rp), INTENT(IN) :: s1,s2
      REAL(rp), INTENT(OUT) :: re(1),se(1)      
      
      INTEGER :: ig
      INTEGER :: it
      INTEGER :: ntry
      INTEGER :: max_it  
      INTEGER :: ndf,dof
      INTEGER :: fail_flag      
      REAL(rp) :: tol   
      REAL(rp) :: dt,t0
      REAL(rp) :: r,s 
      REAL(rp) :: phi(sol%mndof,2),dpdr(sol%mndof,1),dpds(sol%mndof,1)      
      REAL(rp) :: zeta,bathy,xmom,ymom,H
      REAL(rp) :: dzdr,dzds,dbdr,dbds,dHdr,dHds
      REAL(rp) :: dxmdr,dxmds,dymdr,dymds
      REAL(rp) :: uvel,vvel      
      REAL(rp) :: t,f,dfdr,dfds,dfdt      
      
      
      tol = 1d-8
      ntry = 20
      max_it = 1000      
    
      dt = 1d0/real(ntry,rp)
      t0 = 0d0
      fail_flag = 1
guess:DO ig = 1,2*ntry+1
                 
        t = t0                 
        it = 0
 newton:DO 
         
          r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
          s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)
          re(1) = r
          se(1) = s            
                   
          IF (type_flag == 2 .or. type_flag == 4) THEN
                   
            IF (sol%output_type == "adcirc") THEN
              CALL linear_basis(ndf,1,re,se,phi,dpdr,dpds)
            ELSE IF (sol%output_type == "dgswem") THEN
              CALL dgswem_basis(sol%hbp,ndf,1,re,se,phi,dpdr,dpds)   
            ELSE IF (sol%output_type == "dgswe") THEN
              CALL element_basis(et,sol%hbp,ndf,1,re,se,phi,dpdr,dpds)   
            ENDIF 
                                         
            bathy = 0d0
            dbdr = 0d0
            dbds = 0d0
            DO dof = 1,ndf
              bathy = bathy + phi(dof,1)*sol%hb(dof,el,1)
              dbdr = dbdr + dpdr(dof,1)*sol%hb(dof,el,1)
              dbds = dbds + dpds(dof,1)*sol%hb(dof,el,1)
            ENDDO
            IF (type_flag == 2) THEN
              f = bathy - sol_lev
              dfdr = dbdr
              dfds = dbds
            ENDIF
          ENDIF
                   
          IF (type_flag == 3 .or. type_flag == 4) THEN
                   
            IF (sol%output_type == "adcirc") THEN
              CALL linear_basis(ndf,1,re,se,phi,dpdr,dpds)
            ELSE IF (sol%output_type == "dgswem") THEN
              CALL dgswem_basis(sol%p,ndf,1,re,se,phi,dpdr,dpds)    
            ELSE IF (sol%output_type == "dgswe") THEN                
              CALL element_basis(et,sol%p,ndf,1,re,se,phi,dpdr,dpds)    
            ENDIF
                     
            zeta = 0d0
            dzdr = 0d0
            dzds = 0d0
            DO dof = 1,ndf
              zeta = zeta + phi(dof,1)*sol%Z(dof,el,sol_snap)
              dzdr = dzdr + dpdr(dof,1)*sol%Z(dof,el,sol_snap)
              dzds = dzds + dpds(dof,1)*sol%Z(dof,el,sol_snap)
            ENDDO
            IF (type_flag == 3) THEN
              f = zeta - sol_lev
              dfdr = dzdr
              dfds = dzds
            ENDIF
          ENDIF    
                     
          IF (type_flag == 4) THEN
            xmom = 0d0
            ymom = 0d0
            dxmdr = 0d0
            dxmds = 0d0
            dymdr = 0d0
            dymds = 0d0
            DO dof = 1,ndf
              xmom = xmom + phi(dof,1)*sol%Qx(dof,el,sol_snap)
              ymom = ymom + phi(dof,1)*sol%Qy(dof,el,sol_snap)
              dxmdr = dxmdr + dpdr(dof,1)*sol%Qx(dof,el,sol_snap)
              dxmds = dxmds + dpds(dof,1)*sol%Qx(dof,el,sol_snap)
              dymdr = dymdr + dpdr(dof,1)*sol%Qy(dof,el,sol_snap)
              dymds = dymds + dpds(dof,1)*sol%Qy(dof,el,sol_snap)
            ENDDO
            H = zeta + bathy          
            dHdr = dzdr + dbdr
            dHds = dzds + dbds
            uvel = xmom/H
            vvel = ymom/H
            f = uvel**2 + vvel**2 - sol_lev**2
            dfdr = 2d0*(uvel*(dxmdr*H-dHdr*xmom)/H**2 + vvel*(dymdr*H-dHdr*ymom)/H**2)
            dfds = 2d0*(uvel*(dxmds*H-dHds*xmom)/H**2 + vvel*(dymds*H-dHds*ymom)/H**2)
          ENDIF
                   
                   
          dfdt = .5d0*(dfdr*(r2-r1) + dfds*(s2-s1))
                 
          t = t - f/dfdt
          it = it + 1
                   
                   
          IF (abs(f) < tol) THEN


            IF (t <= 1d0+tol .and. t >= -(1d0+tol)) THEN
              r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
              s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)                   
              fail_flag = 0           
!             PRINT*, "Iteration successful"
              EXIT guess
            ELSE                 
!             PRINT("(4(A,I9))"), "  Point outside edge     lev = ", lev, " el = ", el, " tri = ",tri, " vrt = ",vrt
!             PRINT("(A,F12.6,A,I5,A,E24.17,A,E24.17)"), "      t0 = ", t0, "  it = ", it,"  t = ", t, "  f = ", f
!             PRINT*, ""     
              EXIT newton
            ENDIF
          ENDIF
                   
          IF (it > max_it) THEN     
!           PRINT*, "  Max iterations exceeded"
            EXIT newton
          ENDIF
                   
        ENDDO newton
                 
        t0 = t0*(-1d0)                 
        IF (mod(ig,2) == 1) THEN
          t0 = t0 + dt
        ENDIF
                 
      ENDDO guess
                 
      IF (fail_flag == 0) THEN
        re(1) = r
        se(1) = s
      ELSE                
        PRINT*, "  iteration failed ", "r = ",r, " s = ",s 
        re(1) = .5d0*(r1+r2)
        se(1) = .5d0*(s1+s2)
      ENDIF    
      
      
      RETURN
      END SUBROUTINE contour_line_newton      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE contour_line_linear(f1,f2,sol_lev,r1,r2,s1,s2,re,se)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: f1,f2
      REAL(rp), INTENT(IN) :: sol_lev
      REAL(rp), INTENT(IN) :: r1,r2
      REAL(rp), INTENT(IN) :: s1,s2
      REAL(rp), INTENT(OUT) :: re(1),se(1) 
      
      REAL(rp) :: t
      
      t = (2d0*sol_lev-f1-f2)/(f2-f1)
      re(1) = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
      se(1) = .5d0*((1d0-t)*s1 + (1d0+t)*s2)
      
      
      RETURN
      END SUBROUTINE contour_line_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE evaluate_mod
