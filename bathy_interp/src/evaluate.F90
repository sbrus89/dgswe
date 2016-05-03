      MODULE evaluate
      
      IMPLICIT NONE
      
      
      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
 
      SUBROUTINE eval_hb(ele,pte,xe,x,hb)

      USE globals, ONLY: rp,base,eval
      USE find_element, ONLY: in_element
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE transformation, ONLY: element_transformation
      USE bathymetry_interp_mod, ONLY: bathymetry_interp_eval

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ele
      INTEGER, INTENT(IN) :: pte
      REAL(rp), INTENT(IN), OPTIONAL :: xe(2)    
      REAL(rp), INTENT(OUT) :: x(2)
      REAL(rp), INTENT(OUT) :: hb 
      
      INTEGER :: elb,ete,etb,et
      INTEGER :: npt,nd,n
      INTEGER :: info    
      REAL(rp) :: rs(2),r(1),s(1)
      REAL(rp) :: error

      ete = eval%el_type(ele)
      
      IF (PRESENT(xe)) THEN
        ! Take x,y coordinates as evaluation point (ele and pte are unused)
        x(1) = xe(1)      
        x(2) = xe(2)
      ELSE
        ! evaluate x,y coordinates of eval element evaluation point
        CALL element_transformation(eval%nnds(ete),eval%elxy(:,ele,1),eval%elxy(:,ele,2),eval%l(:,pte,ete),x(1),x(2))        
      ENDIF
  
      CALL in_element(x,base%el_type,base%elxy,elb,rs)
          
      etb = base%el_type(elb)  ! element type for base element          
        
      r(1) = rs(1)
      s(1) = rs(2)
      npt = 1  
      CALL shape_functions_area_eval(etb,base%hbp,n,npt,r,s,base%l(:,:,etb))        

               
      ! Evaluate bathymetry at r,s coordinates
      CALL bathymetry_interp_eval(n,base%elhb(:,elb),base%l(:,1,etb),hb)

      RETURN
      END SUBROUTINE eval_hb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      END MODULE evaluate