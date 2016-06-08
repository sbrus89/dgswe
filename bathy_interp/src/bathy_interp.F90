      PROGRAM bathy_interp

      USE globals
      USE allocation
      USE evaluate
      USE basis   
      USE find_element, ONLY: find_element_init
      USE grid_file_mod, ONLY: read_bathy_file
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp
      USE version, ONLY: version_information

      IMPLICIT NONE

      INTEGER :: el_in,el_ex
      INTEGER :: pt_in,pt_ex
      INTEGER :: led_in,led_ex
      INTEGER :: et_in,et_ex
      INTEGER :: nv_in,nv_ex
      INTEGER :: ged,i,nd
      INTEGER :: p,nnds
      INTEGER :: myrank,space
      REAL(rp) :: x(2),xe(2)
      REAL(rp) :: hb

      myrank = 0
      space = 1

      CALL version_information(unit=6) ! print out current git branch/SHA          
      
      CALL read_input()  ! read bathy.inp file      

      CALL sizes(base)    
      CALL sizes(eval)        
      
      CALL read_grid(base)  ! read coarse, fine, and base grids     
      CALL read_grid(eval)
      CALL read_bathy_file(myrank,base%bathy_file,base%hbp,base%ne,base%el_type,nverts,base%depth,base%ect,base%elhb)            
      
      CALL connect(base)
      CALL connect(eval)
      
      mnepn = max(base%mnepn,eval%mnepn)
      
      CALL curvilinear(base)
      CALL curvilinear(eval)       
      
 

      ! evaluate shape functions at evaluation points (to compute r,s -> x,y transformtion)
      CALL shape_functions_eltypes_at_hbp(space,nel_type,eval%np,eval%l,eval%dldr,eval%dlds)
                   

      
           
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      CALL find_element_init(nel_type,nverts,base%np,base%nnds,base%nn,base%xy,base%nepn,base%epn) 

      
      ALLOCATE(base%l(base%mnnds,eval%mnnds,nel_type+2))
      ALLOCATE(eval%elhbxy(eval%mnnds,eval%ne,2))        
      ALLOCATE(eval%hbxy(3,eval%mnnds*eval%ne))          
      

      
      p = eval%hbp
      
      eval%npts = 0
      
      eval%elhb = 0d0
      
      
      PRINT*, " "      
      PRINT*, "Computing vertex points"      
      DO nd = 1,eval%nn
      
        IF (mod(nd,1000) == 0) THEN       
          PRINT*,"Vertex ", nd,"/",eval%nn
        ENDIF        
      
        el_in = eval%epn(1,nd)
        pt_in = 1
        
        xe(1) = eval%xy(1,nd)
        xe(2) = eval%xy(2,nd)
        
        CALL eval_hb(ele=el_in,pte=pt_in,xe=xe,x=x,hb=hb)
        
        eval%npts = eval%npts + 1          
        eval%hbxy(1,eval%npts) = x(1)
        eval%hbxy(2,eval%npts) = x(2)
        eval%hbxy(3,eval%npts) = hb
        
      ENDDO 
      
      
      PRINT*, " "      
      PRINT*, "Computing edge points"
      DO ged = 1,eval%ned 

        IF (mod(ged,1000) == 0) THEN       
          PRINT*,"Edge ", ged,"/",eval%ned
        ENDIF      

        el_in = eval%ged2el(1,ged)
        led_in = eval%ged2led(1,ged)
        et_in = eval%el_type(el_in)
        nv_in = nverts(et_in)
        
        DO i = 1,p-1 
        
          pt_in = mod(led_in,nv_in)*p + i + 1
        
          CALL eval_hb(ele=el_in,pte=pt_in,x=x,hb=hb)
          
          eval%elhb(pt_in,el_in) = hb
          eval%elhbxy(pt_in,el_in,1) = x(1)
          eval%elhbxy(pt_in,el_in,2) = x(2)          
          
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
        nv_in = nverts(et_in)        
        
        IF (mod(et_in,2) == 1) THEN
          nnds = eval%nnds(5)
        ELSE IF (mod(et_in,2) == 0) THEN
          nnds = eval%nnds(6)
        ENDIF
                 
        DO pt_in = nv_in*p+1,nnds
        
          CALL eval_hb(ele=el_in,pte=pt_in,x=x,hb=hb)
          
          eval%elhb(pt_in,el_in) = hb
          eval%elhbxy(pt_in,el_in,1) = x(1)
          eval%elhbxy(pt_in,el_in,2) = x(2)          
          
          eval%npts = eval%npts + 1          
          eval%hbxy(1,eval%npts) = x(1)
          eval%hbxy(2,eval%npts) = x(2)
          eval%hbxy(3,eval%npts) = hb
          
        ENDDO
        
      ENDDO
      
      
      
      ! Populate all vertex values
      DO el_in = 1,eval%ne
        et_in = eval%el_type(el_in)
        nv_in = nverts(et_in)
        
        DO i = 1,nv_in
          nd = (i-1)*p + 1
          pt_in = eval%ect(i,el_in)
          
          eval%elhb(nd,el_in) = eval%hbxy(3,pt_in)  
          eval%elhbxy(nd,el_in,1) = eval%hbxy(1,pt_in) 
          eval%elhbxy(nd,el_in,2) = eval%hbxy(2,pt_in)                
          
        ENDDO
      ENDDO
      
      
      
      
      ! Populate all edge values      
      DO ged = 1,eval%ned            
      
        IF (eval%ed_type(ged) == 0) THEN
        
          el_in = eval%ged2el(1,ged)
          led_in = eval%ged2led(1,ged)
          et_in = eval%el_type(el_in)
          nv_in = nverts(et_in)        
                      
          el_ex = eval%ged2el(2,ged)
          led_ex = eval%ged2led(2,ged)
          et_ex = eval%el_type(el_ex)
          nv_ex = nverts(et_ex)
          
          DO i = 1,p-1
          
            pt_in = mod(led_in,nv_in)*p + i + 1
            pt_ex = mod(led_ex,nv_ex)*p + p - i + 1                      

            
            eval%elhb(pt_ex,el_ex) = eval%elhb(pt_in,el_in)
            eval%elhbxy(pt_ex,el_ex,1) = eval%elhbxy(pt_in,el_in,1) 
            eval%elhbxy(pt_ex,el_ex,2) = eval%elhbxy(pt_in,el_in,2)              
            
          ENDDO
          
        ENDIF      
        
      ENDDO
      
      CALL write_results(eval,base)
      
      CALL interp_check()          
          
      END PROGRAM bathy_interp