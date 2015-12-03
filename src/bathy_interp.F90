      PROGRAM bathy_interp

      USE globals
      USE allocation
      USE read_grid
      USE evaluate
      USE basis   
      USE find_element, ONLY: in_element

      IMPLICIT NONE

      INTEGER :: ele,elb
      INTEGER :: ete,etb,et
      INTEGER :: npts,pte,nd,n
      INTEGER :: info
      REAL(rp) :: x(2)
      REAL(rp) :: xe(1),ye(1)
      REAL(rp) :: r(1),s(1)
      REAL(rp) :: hb
 

      CALL version () ! print out current git branch/SHA
      
      CALL read_input()  ! read bathy.inp file      
      
      CALL read_grids()  ! read coarse, fine, and base grids     
      
      CALL connect(base)
      CALL connect(eval)
      
      mnepn = max(base%mnepn,eval%mnepn)
      
      CALL vandermonde(base) ! compute Vandermonde matricies
      CALL vandermonde(eval)
            
      CALL re_vert(base) ! find reference element verticies
      CALL re_vert(eval)     
      
      CALL curvilinear(base)
      CALL curvilinear(eval)       
      
      CALL eval_pts()      ! get area evaluation points   
      CALL function_eval() ! evaluate basis and shape functions for eval grid evaluation points            
                   
      CALL read_bathy()  
      
           
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      tree_xy  => kdtree2_create(base%vxy, rearrange=.true., sort=.true.)  
      ALLOCATE(closest(srchdp))   

      
      ALLOCATE(base%l(base%mnnds,mnept,nel_type+2))

      DO ele = 1,eval%ne  ! Calculate error integral in the fine grid elements

        ete = eval%el_type(ele) ! element type for evel element
        npts = nept(ete)        ! number of evaluation points
        
        DO pte = 1,npts                       
        
          ! evaluate x,y coordinates of fine element quadrature points
          x(1) = 0d0
          x(2) = 0d0

          DO nd = 1,eval%nnds(ete)
            x(1) = x(1) + eval%l(nd,pte,ete)*eval%elxy(nd,ele,1)  
            x(2) = x(2) + eval%l(nd,pte,ete)*eval%elxy(nd,ele,2)
          ENDDO   
          
          PRINT*,ele,pte
          PRINT*,x(1),x(2)
          CALL in_element(x(1:2),elb)
          
          etb = base%el_type(elb)  ! element type for base element          
        
          ! calculate r,s coordinates of fine element quadrature points
          CALL newton(base,x(1),x(2),1,elb,r,s)  
               
          IF (mod(etb,2) == 1) THEN
            et = 5
            CALL tri_basis(base%hbp,n,1,r,s,base%l(:,:,et))
          ELSE IF (mod(etb,2) == 0) THEN
            et = 6
            CALL quad_basis(base%hbp,n,1,r,s,base%l(:,:,et))
          ENDIF       
        
          CALL DGETRS("N",n,1,base%V(1,1,et),base%mnnds,base%ipiv(1,et),base%l(1,1,et),base%mnnds,info)                          
               
               
          ! Evaluate bathymetry at r,s coordinates
          hb = 0d0
          DO nd = 1,n
            hb = hb + base%l(nd,1,et)*base%elhb(nd,elb)
          ENDDO            

          eval%elhb(pte,ele) = hb
          
          
        ENDDO 
        
      ENDDO 
      
          
      END PROGRAM bathy_interp