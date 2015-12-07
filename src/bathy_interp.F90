      PROGRAM bathy_interp

      USE globals
      USE allocation
      USE read_grid
      USE evaluate
      USE basis   
      USE find_element, ONLY: in_element

      IMPLICIT NONE

      INTEGER :: el_in,el_ex
      INTEGER :: pt_in,pt_ex
      INTEGER :: led_in,led_ex
      INTEGER :: et_in,et_ex
      INTEGER :: nv_in,nv_ex
      INTEGER :: ged,i
      INTEGER :: p,nnds
      REAL(rp) :: x(2)
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
      
      

      
      p = eval%hbp
      
      eval%npts = 0
      
      PRINT*, "Computing edge and vertex points"
      DO ged = 1,eval%ned 

        IF (mod(ged,1000) == 0) THEN       
          PRINT*,"Edge ", ged,"/",eval%ned
        ENDIF      

        el_in = eval%ged2el(1,ged)
        led_in = eval%ged2led(1,ged)
        et_in = eval%el_type(el_in)
        nv_in = eval%nverts(et_in)
        
        DO i = 1,p 
        
          pt_in = mod(led_in,nv_in)*p + i
        
          CALL eval_hb(el_in,pt_in,x,hb)
          
          eval%elhb(pt_in,el_in) = hb

          eval%npts = eval%npts + 1          
          eval%hbxy(1,eval%npts) = x(1)
          eval%hbxy(2,eval%npts) = x(2)
          eval%hbxy(3,eval%npts) = hb
          

            
        ENDDO 
        
      ENDDO 
      
      PRINT*, " "
      PRINT*, "Computing interior points"      
      DO el_in = 1,eval%ne
      
        IF (mod(el_in,1000) == 0) THEN       
          PRINT*,"Element ", el_in,"/",eval%ne
        ENDIF            
      
        et_in = eval%el_type(el_in)
        nv_in = eval%nverts(et_in)        
        
        IF (mod(et_in,2) == 1) THEN
          nnds = eval%nnds(5)
        ELSE IF (mod(et_in,2) == 0) THEN
          nnds = eval%nnds(6)
        ENDIF
                 
        DO pt_in = nv_in*p+1,nnds
        
          CALL eval_hb(el_in,pt_in,x,hb)
          
          eval%elhb(pt_in,el_in) = hb
          
          eval%npts = eval%npts + 1          
          eval%hbxy(1,eval%npts) = x(1)
          eval%hbxy(2,eval%npts) = x(2)
          eval%hbxy(3,eval%npts) = hb
          
        ENDDO
        
      ENDDO
      
      
      DO ged = 1,eval%ned            
      
        IF (eval%bed_flag(ged) == 0) THEN
        
          el_in = eval%ged2el(1,ged)
          led_in = eval%ged2led(1,ged)
          et_in = eval%el_type(el_in)
          nv_in = eval%nverts(et_in)        
                      
          el_ex = eval%ged2el(2,ged)
          led_ex = eval%ged2led(2,ged)
          et_ex = eval%el_type(el_ex)
          nv_ex = eval%nverts(et_ex)
          
          DO i = 1,p+1
          
            pt_in = mod(led_in,nv_in)*p + i
            pt_ex = mod(led_ex,nv_ex)*p + p - i + 2
            
            IF (pt_in == nv_in*p) THEN
              pt_in = 1
            ENDIF
            
            IF (pt_ex == nv_ex*p) THEN
              pt_ex = 1
            ENDIF
            
            eval%elhb(pt_ex,el_ex) = eval%elhb(pt_in,el_in)
            
          ENDDO
          
        ENDIF      
        
      ENDDO
      
      CALL write_results(eval)
      
          
      END PROGRAM bathy_interp