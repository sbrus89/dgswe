      MODULE evaluate_mod

      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type

      IMPLICIT NONE


      CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_solution(npt,r,s,fig)
      
      IMPLICIT NONE
      

      

      
      fig%snap_min = 1d10
      fig%snap_max = -1d10
      
      
          IF (fig%type_flag == 2 .or. fig%type_flag == 4) THEN
          
            CALL element_basis(et,p,ndf,npplt(i),r(:,i),s(:,i),fig%phi(:,:,et))            
          
            DO nd = 1,npplt(i)
              hb_val(nd) = 0d0
              DO dof = 1,ndf
                hb_val(nd) = hb_val(nd) + fig%phi(dof,nd,et)*hb(dof,el,1)
              ENDDO
            ENDDO
            
            IF (fig%type_flag == 2) THEN
              DO nd = 1,npplt(i) 
                fig%sol_val(nd,el) = hb_val(nd)
              ENDDO
            ENDIF
          ENDIF                         
     
          IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN
          
            CALL element_basis(et,p,fig%ndof(et),npplt(i),r(:,i),s(:,i),fig%phi(:,:,et))           
          
            DO nd = 1,npplt(i)
              zeta_val(nd) = 0d0
              DO dof = 1,fig%ndof(et)
                zeta_val(nd) = zeta_val(nd) + fig%phi(dof,nd,et)*Z(dof,el,snap)
              ENDDO
            ENDDO
            
            IF (fig%type_flag == 3) THEN
              DO nd = 1,npplt(i)
                fig%sol_val(nd,el) = zeta_val(nd)
              ENDDO
            ENDIF
          ENDIF
          
          IF (fig%type_flag == 4) THEN
            DO nd = 1,npplt(i)
              Qx_val(nd) = 0d0
              Qy_val(nd) = 0d0
              DO dof = 1,fig%ndof(et)
                Qx_val(nd) = Qx_val(nd) + fig%phi(dof,nd,et)*Qx(dof,el,snap)
                Qy_val(nd) = Qy_val(nd) + fig%phi(dof,nd,et)*Qy(dof,el,snap)                
              ENDDO
            ENDDO
            
            DO nd = 1,npplt(i) 
              H = zeta_val(nd)+hb_val(nd)
              u_vel = Qx_val(nd)/H
              v_vel = Qy_val(nd)/H
              fig%sol_val(nd,el) = sqrt(u_vel**2+v_vel**2)
            ENDDO
          ENDIF                

  
      
!           IF (fig%type_flag == 2 .and. fig%sol_val(nd,el) < fig%h0) THEN
!             PRINT "(A,I7,A,F15.7)", "Element depth is less than tolerance: ", el, " hb = ", fig%sol_val(nd,el)          
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