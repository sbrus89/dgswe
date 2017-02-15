      PROGRAM spline

      USE globals, ONLY: rp,base,eval,nverts,nel_type, &
                         ax,bx,cx,dx,ay,by,cy,dy,dt, &
                         theta_tol,sig, &
                         nfbnds,fbnds,fbnds_xy,nfbnd2el,fbnd2el
      USE allocation, ONLY: sizes
                         
      USE calc_spline, ONLY: calc_cubic_spline,eval_cubic_spline, &
                             newton,spline_init,evaluate,update_elxy_spline
      USE check, ONLY: check_angle,check_deformation,check_transformations, &
                       l2_project,quad_interp
      USE find_element, ONLY: in_element,check_elem,find_element_init,sub_element 
      USE curvilinear_nodes_mod, ONLY: shape_functions_linear_at_ctp,eval_coordinates_curved
      USE version, ONLY: version_information

      IMPLICIT NONE
      INTEGER :: i,j,k,n,bou,num,nmax,qpts,ebou
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert,nv
      INTEGER :: el_in,found
      INTEGER :: btype
      INTEGER :: n1,n2
      INTEGER :: base_bed,base_bou,base_led,base_vert
      INTEGER :: neval,nbase
      INTEGER :: error_flag,ntry,try,success
      INTEGER, ALLOCATABLE, DIMENSION(:) :: base_els
      REAL(rp) :: htest,ti,tpt,r,rpt,ra,xs,ys,r0,dr
      REAL(rp) :: d1,d2,d3,t1,t2,xr(2),xa(2),rs(2)
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: theta1,theta2
      REAL(rp) :: diff
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: x,y
      INTEGER :: eds(4)
      INTEGER :: all_bnd
      REAL(rp), PARAMETER :: it_tol = 1d-3

      


      PRINT "(A)", " "
      
      CALL version_information(unit=6)
      
      CALL read_input()
      
      CALL sizes(base)
      CALL sizes(eval)

      CALL read_grid(base)
      CALL read_grid(eval)
      
      CALL connect(base)
      CALL connect(eval)
      
      sig = 0d0
      
      CALL shape_functions_at_qpts(base)
      CALL shape_functions_linear_at_ctp(nel_type,base%np,base%psiv)      
      
      CALL spline_init(num,nmax)
      
      CALL find_element_init(nel_type,nverts,eval%np,eval%nnds,nfbnds,fbnds_xy,nfbnd2el,fbnd2el)
      


    
      OPEN(unit=30,file='spline.out')   
      OPEN(unit=90,file='max_deform.out')
    
      WRITE(30,*) num
      WRITE(90,*) num
      
 

 calc:DO bou = 1,base%nbou
      
        btype = base%fbseg(2,bou)               
        
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
        
          n = base%fbseg(1,bou)    ! n nodes, n-1 subintervals
          
          PRINT "(A)", " "          
          PRINT "(A)", " "          
          PRINT "(A)", " "          
          PRINT "(A)", " "
          PRINT "(A,I5)", "Normal flow boundary ",bou
          PRINT "(A,I5)", "Normal flow boundary nodes ",n
 


          !!!!!!!!!!!!!!!!!!!
          ! x value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(1,bou,n,btype,sig,ax(:,bou),bx(:,bou), &
                                                   cx(:,bou),dx(:,bou), &
                                                   dt(:,bou))
       

          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(2,bou,n,btype,sig,ay(:,bou),by(:,bou), &
                                                   cy(:,bou),dy(:,bou), &
                                                   dt(:,bou))
          
                  
           
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! modify spline if necessary
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WRITE(90,"(I6)") n-1
          
          DO nd = 1,n-1
            PRINT*, "node: ",nd ,"node #:", base%fbnds(nd,bou)
          
            ti = 0d0        ! find starting parameter value for found edge
            DO j = 1,nd-1
              ti = ti + dt(j,bou)
            ENDDO            
                                            
            CALL check_angle(bou,nd,dt(nd,bou),ti,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                  ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou), &
                                                  theta1,theta2)            
                                              
            CALL check_deformation(bou,nd,dt(nd,bou),ti,theta1,theta2,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                                      ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou)) 
                                                                      
            CALL update_elxy_spline(base,nverts,bou,nd,dt(nd,bou),ti,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                                       ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou))                                                                                  
                                                        
            PRINT*, "------------------"       
            
!             IF(base%fbnds(nd,bou) == 1314) THEN
!               STOP
!             ENDIF
                        
          ENDDO          



          CALL write_spline(n,bou,base%fbnds(:,bou))   
          
!           PAUSE
          
        ENDIF
      ENDDO calc
        
        
      CALL check_transformations(base,nverts)  
        
        
      DO i = 1,10  
        PRINT*,""  
      ENDDO
      
!       PAUSE
      
      
      
      
      ALLOCATE(x(eval%ctp+1),y(eval%ctp+1)) 
      ALLOCATE(eval%bndxy(2,eval%ctp+1,eval%nvel,eval%nbou))
      
      CALL shape_functions_at_qpts(eval)
      CALL shape_functions_linear_at_ctp(nel_type,eval%np,eval%psiv)          
      
      
     
      
      OPEN(unit=60,file='eval_nodes.out')       
      WRITE(60,*) eval%nbou      
          
      DO ebou = 1,eval%nbou
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! evaluate spline at eval grid points
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        
        neval = eval%fbseg(1,ebou)
        btype = eval%fbseg(2,ebou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
             btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries          
          
          
          WRITE(60,*) eval%ctp*(neval-1) + 1

          DO i = 1,neval-1  
             
            n1 = eval%fbnds(i,ebou)
            n2 = eval%fbnds(i+1,ebou)
                   
            n1x = eval%xy(1,n1)
            n1y = eval%xy(2,n1)
          
            n2x = eval%xy(1,n2)
            n2y = eval%xy(2,n2)             
            
            IF (i == neval-1) THEN
              n = eval%ctp
            ELSE 
              n = eval%ctp-1
            ENDIF
                   
!             PRINT*, "EVALUATING SPLINE COORDINATES"        
            DO j = 0,n                                     
            
              rpt = eval%rpts(j+1)   
              
              IF (j == 0) THEN
                ra = rpt + 1d-2          ! add an offset avoid ambiguity with verticies
              ELSE IF (j == eval%ctp) THEN
                ra = rpt - 1d-2
              ELSE 
                ra = rpt
              ENDIF
                                                            
              xa(1) = .5d0*(1d0-ra)*n1x + .5d0*(1d0+ra)*n2x 
              xa(2) = .5d0*(1d0-ra)*n1y + .5d0*(1d0+ra)*n2y                                  
              
              PRINT*, "FINDING ELEMENT FOR POINT: ", i, " NODE: ",n1
              success = 0
              CALL in_element(xa,base%el_type,base%elxy,el_in,rs,eds,base_vert,base_els) 
              
    try_elems:DO el = 1,size(base_els) ! loop over elements in case a point isn't found in the -1,1 interval for the closest element by area difference
    
                el_in = base_els(el)
                PRINT*, "USING ELEMENT ", el_in
                CALL sub_element(xa,el_in,base%el_type,base%elxy,diff,eds,base_vert,rs) ! this is a little sloppy since it's a reapeat call for the first element
                                                                                        ! but it's needed to update information for elements in the loop                                                                                        
                CALL find_edge(n1,n2,xa,el_in,eds,base_bou,base_bed,base_led)  ! find base edge (to get correct spline coefficients) 
              
              
                ! take care of case where a channel spanned by a single triangluar element           
                ! meaning all three element nodes that are on spline boundaries             
                ! this can cause an element to be chosen when its boundary edge is on the other side of the channel     
                ! particularly when ctp=2 and the eval grid is a derefinement of the base grid  
                
                IF (base%el_type(el_in) == 1) THEN                   
                  all_bnd = 0                                       
          search: DO k = 1,nfbnds                                    
                    IF (base%ect(base_led,el_in) == fbnds(k)) THEN   ! check to see if node that is not on the found edge is also on a spline boundary
                        all_bnd = 1
                       EXIT search            
                    ENDIF
                  ENDDO  search                 
                    
                  IF (all_bnd == 1) THEN                               
                    IF (base_led == eds(3)) THEN                     ! check if this edge is the farthest from the the xa point
                      CYCLE try_elems                                ! if it is, this element should be skipped 
                    ENDIF                                            ! because it has a boundary edge on the other side of the channel
                  ENDIF
                ENDIF
              
                nd = base_bed 
                bou = base_bou
              
                xr(1) = .5d0*(1d0-rpt)*n1x + .5d0*(1d0+rpt)*n2x
                xr(2) = .5d0*(1d0-rpt)*n1y + .5d0*(1d0+rpt)*n2y                
              
! !                 check to make sure vertex offset used to find element isn't 
! !                 too large that the found element isn't connected to the vertex point
!                 IF (j == 0 .OR. j == eval%ctp) THEN           
!                   CALL check_elem(xr,el_in)         
!                 ENDIF               
            
                ti = 0d0        ! find starting parameter value for found edge
                DO k = 1,nd-1
                  ti = ti + dt(k,bou)
                ENDDO              
                      
                   
                IF (base_vert /= 0) THEN  ! if the eval point is near a base vertex, use start of segement              
                  r0 = -1d0               ! as the initial guess
                ELSE
                  SELECT CASE (base_led)  ! use the initial guess that corresponds to the closest edge
                    CASE(1)         
                      r0 = rs(1)
                    CASE(2) 
                      r0 = rs(2)
                    CASE(3) 
                      r0 = rs(1)
                  END SELECT
                ENDIF
              
              
                   
                r = r0     
                CALL evaluate(r,dt(nd,bou),ti,xr,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                 ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou), &
                                                 x(j+1),y(j+1),error_flag)

!                 WRITE(60,*) x(j+1),y(j+1)                                                  
              
                           
                IF (abs(r)-1d0 > it_tol) THEN
                  PRINT "(A,F24.17)", "ERROR: R VALUE NOT FOUND IN INTERVAL, R = ", r
                  PRINT "(2(F24.17))", abs(r)-1d0, it_tol
!                   STOP
                ELSE 
                  PRINT "(A,F24.17)", "R = ", r 
                  success = 1
                  EXIT try_elems
                ENDIF            
                

              
              ENDDO try_elems
              
              IF (success == 1) THEN
                WRITE(60,*) x(j+1),y(j+1)              
              ELSE 
                PRINT "(A,F24.17)", "ERROR: R VALUE NOT FOUND IN INTERVAL FOR ANY NEIGHBOR ELEMENTS"
                STOP
              ENDIF
              
          

                
              
              IF (j == 0) THEN
                eval%xy(1,n1) = x(j+1)
                eval%xy(2,n1) = y(j+1)
              ELSE IF (j == eval%ctp) THEN  
                eval%xy(1,n2) = x(j+1)
                eval%xy(2,n2) = y(j+1)            
              ELSE 
                eval%bndxy(1,j,i,ebou) = x(j+1)
                eval%bndxy(2,j,i,ebou) = y(j+1)                
              ENDIF              
              
                       
            ENDDO  
            
!               IF (i == 17) THEN
!                 STOP
!               ENDIF            

!                 IF (n1 == 72326) THEN         
!                   STOP
!                 ENDIF    
            
            WRITE(40,"(I8,1X,10(E24.17,1X))") n1, (x(j), j=1,eval%ctp)
            WRITE(40,"(I8,1X,10(E24.17,1X))") n1, (y(j), j=1,eval%ctp)
            

            IF (i == neval-1) THEN
              WRITE(40,"(I8,1X,10(E24.17,1X))") n2, x(eval%ctp+1)
              WRITE(40,"(I8,1X,10(E24.17,1X))") n2, y(eval%ctp+1)            
            ENDIF
            
            
            PRINT*, "------------------------------------------------------------"            
            PRINT*, " "
            
            
          ENDDO

        ELSE
        

          WRITE(60,*) 0
              
        ENDIF 
        

      ENDDO
      CLOSE(30)
      CLOSE(60)
      
      
      PRINT*, "Writing .cb file"
      CALL write_cb_file(eval)          
      
      PRINT*, "Checking coordniate transformations"
      CALL eval_coordinates_curved(eval%ctp,eval%nnds,nverts,eval%el_type_spline,eval%xy,eval%ect,eval%fbseg,eval%fbnds, &
                                   eval%nnfbed,eval%nfbedn,eval%nfbednn,eval%ged2el,eval%ged2led, &
                                   eval%psiv,eval%bndxy,eval%elxy_spline) 
                                              
      
      CALL check_transformations(eval,nverts)  
      
      
      PRINT*, "Writing grid output file"
      CALL write_grid(eval,base%grid_file)           
     

      END PROGRAM spline
      
      

      
      
