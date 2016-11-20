      MODULE evaluate_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type

      IMPLICIT NONE


      CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_solution(et,sol_type,npts,r,s,snap,sol_val)

      USE globals, ONLY: hb,Z,Qx,Qy
      USE read_dginp, ONLY: p,hbp
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: sol_type
      INTEGER, INTENT(IN) :: npts
      REAL(rp), DIMENSION(:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:), INTENT(IN) :: s
      INTEGER, INTENT(IN) :: snap
      REAL(rp), DIMENSION(:), INTENT(OUT) :: sol_val

      INTEGER :: pt,dof
      INTEGER :: mndf,ndf
      REAL(rp) :: H,u_vel,v_val
      REAL(rp) :: hb_val(npts),zeta_val(npts)
      REAL(rp) :: Qx_val(npts),Qy_val(npts)
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: phi

      mndf = (p+1)**2
      ALLOCATE(phi(mndf,npts))
      
      
      IF (sol_type == 2 .or. sol_type == 4) THEN
          
        CALL element_basis(et,hbp,ndf,npts,r,s,phi)            
          
        DO pt = 1,npts
          hb_val(pt) = 0d0
          DO dof = 1,ndf
            hb_val(pt) = hb_val(pt) + phi(dof,pt)*hb(dof,el,1)
          ENDDO
        ENDDO
            
        IF (sol_type == 2) THEN
          DO pt = 1,npts 
            sol_val(pt) = hb_val(pt)
          ENDDO
        ENDIF
      ENDIF                         
     
      IF (sol_type == 3 .or. sol_type == 4) THEN
          
        CALL element_basis(et,p,ndf,npts,r,s,phi)           
          
        DO pt = 1,npts
          zeta_val(pt) = 0d0
          DO dof = 1,ndf
            zeta_val(pt) = zeta_val(pt) + phi(dof,pt)*Z(dof,el,snap)
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
          DO dof = 1,ndf
            Qx_val(pt) = Qx_val(pt) + phi(dof,pt)*Qx(dof,el,snap)
            Qy_val(pt) = Qy_val(pt) + phi(dof,pt)*Qy(dof,el,snap)                
          ENDDO
        ENDDO
          
        DO pt = 1,npts 
          H = zeta_val(pt)+hb_val(pt)
          u_vel = Qx_val(pt)/H
          v_vel = Qy_val(pt)/H
          sol_val(pt) = sqrt(u_vel**2+v_vel**2)
        ENDDO
      ENDIF                

  
      
!           IF (sol_type == 2 .and. sol%sol_val(pt,el) < sol%h0) THEN
!             PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb = ", sol%sol_val(pt,el)          
!           ENDIF           
      
      RETURN
      END SUBROUTINE evaluate_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE evaluate_basis(mnpp,mndof,nel_type,npplt,r,s,sol)
      
      USE basis, ONLY: element_nodes,element_basis,linear_basis          
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: mnpp
      INTEGER, INTENT(IN) :: mndof
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: s
      TYPE(plot_type), INTENT(INOUT) :: sol
      
      INTEGER :: i,j
      INTEGER :: et
      INTEGER :: n

      
      ALLOCATE(sol%phi(mndof,mnpp,nel_type))    

      DO et = 1,nel_type

#ifndef adcirc      
        CALL element_basis(et,sol%p,sol%ndof(et),npplt(et),r(:,et),s(:,et),sol%phi(:,:,et))       
#else        
        CALL linear_basis(npplt(et),r(:,et),s(:,et),sol%phi(:,:,et))
        sol%ndof(et) = 3
#endif
        
      ENDDO           
      
      RETURN
      END SUBROUTINE evaluate_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE find_solution_minmax(ne,el_type,el_in,npplt,fig,Z,hb,Qx,Qy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt      
      TYPE(plot_type), INTENT(INOUT) :: fig
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: Z           
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: hb      
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Qx     
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Qy      
      
      REAL(rp) :: temp_min,temp_max
      
      
      IF (fig%type_flag == 3) THEN       
        CALL evaluate_depth_solution(ne,el_type,el_in,npplt,Z,fig)      
      ENDIF
      
      IF (fig%plot_sol_option /= 1) THEN
        RETURN
      ENDIF

      IF (fig%type_flag == 4 .AND. fig%cscale_option == "auto-all") THEN      
        CALL evaluate_velocity_solution(ne,el_type,el_in,npplt,Qx,Qy,Z,hb,fig)
      ENDIF      
      
      
      IF (fig%snap_min < fig%sol_min) THEN
        fig%sol_min = fig%snap_min
       ENDIF
          
      IF (fig%snap_max > fig%sol_max) THEN
        fig%sol_max = fig%snap_max
      ENDIF      
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



      END MODULE evaluate_mod