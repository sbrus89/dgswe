      MODULE evaluate_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type

      IMPLICIT NONE


      CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_depth_solution(ne,el_type,el_in,npplt,H,fig)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: H
      TYPE(plot_type), INTENT(INOUT) :: fig
      
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(fig%sol_val)) THEN
        mnpp = MAXVAL(npplt)
        ALLOCATE(fig%sol_val(mnpp,ne)) 
      ENDIF
      
      fig%snap_min = 1d10
      fig%snap_max = -1d10

 elem:DO el = 1,ne
 
        IF(el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        npts = npplt(et)
        ndf = fig%ndof(et)

        
        DO nd = 1,npts
          fig%sol_val(nd,el) = 0d0
          DO dof = 1,ndf
            fig%sol_val(nd,el) = fig%sol_val(nd,el) + H(dof,el)*fig%phi(dof,nd,et)           
          ENDDO    
                  
          IF (fig%sol_val(nd,el) < fig%snap_min) THEN
            fig%snap_min = fig%sol_val(nd,el)
          ENDIF
          
          IF (fig%sol_val(nd,el) > fig%snap_max) THEN
            fig%snap_max = fig%sol_val(nd,el)
          ENDIF
          
          IF (fig%type_flag == 2 .and. fig%sol_val(nd,el) < fig%h0) THEN
            PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb = ", fig%sol_val(nd,el)          
          ENDIF          
        ENDDO       
        
      ENDDO elem      
      
      RETURN
      END SUBROUTINE evaluate_depth_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate_velocity_solution(ne,el_type,el_in,npplt,Qx,Qy,Z_val,hb_val,vel)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: Qx
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: Qy      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: Z_val
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: hb_val  
      TYPE(plot_type), INTENT(INOUT) :: vel      
      
      REAL(rp) :: Qx_val,Qy_val,H_val
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(vel%sol_val)) THEN
        mnpp = MAXVAL(npplt)
        ALLOCATE(vel%sol_val(mnpp,ne)) 
      ENDIF
      
      vel%snap_min = 1d10
      vel%snap_max = -1d10

 elem:DO el = 1,ne
 
        IF(el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        npts = npplt(et)
        ndf = vel%ndof(et)
        
        DO nd = 1,npts
          Qx_val = 0d0
          Qy_val = 0d0
          DO dof = 1,ndf
            Qx_val = Qx_val + Qx(dof,el)*vel%phi(dof,nd,et)
            Qy_val = Qy_val + Qy(dof,el)*vel%phi(dof,nd,et)            
          ENDDO            
          H_val = Z_val(nd,el) + hb_val(nd,el)
          vel%sol_val(nd,el) = sqrt((Qx_val/H_val)**2 + (Qy_val/H_val)**2)   
          
          IF(vel%sol_val(nd,el) < vel%snap_min) THEN
            vel%snap_min = vel%sol_val(nd,el)
          ENDIF
          
          IF (vel%sol_val(nd,el) > vel%snap_max) THEN
            vel%snap_max = vel%sol_val(nd,el)
          ENDIF
        ENDDO
        
      ENDDO elem      
      
      
      RETURN
      END SUBROUTINE evaluate_velocity_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate_basis(p,mnpp,mndof,nel_type,npplt,r,s,sol)
      
      USE basis, ONLY: element_nodes,element_basis,linear_basis          
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p
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
        CALL element_basis(et,p,sol%ndof(et),npplt(et),r(:,et),s(:,et),sol%phi(:,:,et))       
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