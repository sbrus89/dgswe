      PROGRAM spline

      USE globals, ONLY: rp,base,eval,ctp,nverts,np,nnds,nel_type, &
                         ax,bx,cx,dx,ay,by,cy,dy,dt, &
                         rpts,theta_tol,sig, &
                         nfbnds,fbnds_xy,nfbnd2el,fbnd2el
      USE allocation, ONLY: sizes
                         
      USE calc_spline, ONLY: calc_cubic_spline,eval_cubic_spline, &
                             newton,spline_init
      USE check, ONLY: check_angle,check_deformation,l2_project,quad_interp
      USE find_element, ONLY: in_element,check_elem,find_element_init 

      IMPLICIT NONE
      INTEGER :: i,j,k,n,bou,num,nmax,qpts,ebou
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert
      INTEGER :: el_in,found
      INTEGER :: btype
      INTEGER :: n1,n2
      INTEGER :: base_bed,base_bou
      INTEGER :: neval,nbase
      REAL(rp) :: htest,ti,tpt,r,ra,xs,ys
      REAL(rp) :: d1,d2,d3,t1,t2,xr(2),xa(2)
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: theta1,theta2
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: x,y
      INTEGER :: eind
      INTEGER :: eds(4)
      CHARACTER(100) :: name
      CHARACTER(1) :: ctp_char

      


      PRINT "(A)", " "
      
      CALL version()
      
      CALL read_input()
      
      CALL sizes()

      CALL read_grid(base)
      CALL read_grid(eval)
      
      CALL connect(base)
      CALL connect(eval)
      
      sig = 0d0
      
      
      CALL spline_init(num,nmax)
      
      CALL find_element_init(nel_type,nverts,np,nnds,nfbnds,fbnds_xy,nfbnd2el,fbnd2el)
      
      
      ALLOCATE(x(ctp+1),y(ctp+1))

    
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


!           PAUSE
 


          !!!!!!!!!!!!!!!!!!!
          ! x value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(1,bou,n,sig,ax(:,bou),bx(:,bou), &
                                             cx(:,bou),dx(:,bou), &
                                             dt(:,bou))


          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(2,bou,n,sig,ay(:,bou),by(:,bou), &
                                             cy(:,bou),dy(:,bou), &
                                             dt(:,bou))
          
                  
           
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! modify spline if necessary
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WRITE(90,*) n-1
          
          DO nd = 1,n-1
          
            ti = 0d0        ! find starting parameter value for found edge
            DO j = 1,nd-1
              ti = ti + dt(j,bou)
            ENDDO            
                                            
            CALL check_angle(bou,nd,dt(nd,bou),ti,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                  ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou))            
                                              
            CALL check_deformation(bou,nd,dt(nd,bou),ti,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                                        ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou))                           
                        
          ENDDO          



          CALL write_spline(n,bou)      
          
        ENDIF
      ENDDO calc
        
      DO i = 1,10  
        PRINT*,""  
      ENDDO
      
      
      
      
      
      
      eind = INDEX(ADJUSTL(TRIM(eval%grid_file)),".",.false.)   
      name = ADJUSTL(TRIM(eval%grid_file(1:eind-1)))
      WRITE(ctp_char,"(I1)") ctp      
      
      OPEN(unit=60,file='eval_nodes.out')       
      OPEN(unit=40,file=ADJUSTL(TRIM(name)) // "_ctp" // ctp_char // ".cb")    
      
      WRITE(40,"(I8,19x,A)") eval%nbou, "! total number of normal flow boundaries"
      WRITE(40,"(2(I8),19x,A)") eval%nvel,ctp, "! max number of normal flow nodes, ctp order"  
      
      WRITE(60,*) eval%nbou      
          
      DO ebou = 1,eval%nbou
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! evaluate spline at eval grid points
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
        
        neval = eval%fbseg(1,ebou)
        btype = eval%fbseg(2,ebou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
             btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries          
          
          WRITE(40,"(2(I8),19x,A)") neval,btype, "! number of nodes in boundary, boundary type"            
          WRITE(60,*) ctp*(neval-1) + 1

          DO i = 1,neval-1  
             
            n1 = eval%fbnds(i,ebou)
            n2 = eval%fbnds(i+1,ebou)
                   
            n1x = eval%xy(1,n1)
            n1y = eval%xy(2,n1)
          
            n2x = eval%xy(1,n2)
            n2y = eval%xy(2,n2)             
            
            IF (i == neval-1) THEN
              n = ctp
            ELSE 
              n = ctp-1
            ENDIF
                   
!             PRINT*, "EVALUATING SPLINE COORDINATES"        
            DO j = 0,n                                     
            
              r = rpts(j+1)   
              
              IF (j == 0) THEN
                ra = r + 1d-2          ! add an offset avoid ambiguity with verticies
              ELSE IF (j == ctp) THEN
                ra = r - 1d-2
              ELSE 
                ra = r
              ENDIF
                                                            
              xa(1) = .5d0*(1d0-ra)*n1x + .5d0*(1d0+ra)*n2x 
              xa(2) = .5d0*(1d0-ra)*n1y + .5d0*(1d0+ra)*n2y                                  
              
              PRINT*, "FINDING ELEMENT FOR POINT: ", i, " NODE: ",n1
              CALL in_element(xa,base%el_type,base%elxy,el_in,eds)   
              CALL find_edge(n1,n2,xa,el_in,eds,base_bou,base_bed)  ! find base edge (to get correct spline coefficients) 
              
              nd = base_bed 
              bou = base_bou
              
              xr(1) = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
              xr(2) = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y                
              
              ! check to make sure vertex offset used to find element isn't 
              ! too large that the found element isn't connected to the vertex point
              IF (j == 0 .OR. j == ctp) THEN           
                CALL check_elem(xr,el_in)         
              ENDIF               
            
              ti = 0d0        ! find starting parameter value for found edge
              DO k = 1,nd-1
                ti = ti + dt(k,bou)
              ENDDO              
                      
                   
              CALL newton(r,dt(nd,bou),ti,xr,ax(nd,bou),bx(nd,bou),cx(nd,bou),dx(nd,bou), &
                                             ay(nd,bou),by(nd,bou),cy(nd,bou),dy(nd,bou), &
                                             x(j+1),y(j+1))
              

!               ! Try new initial guess if minimum was not found in (-1,1) interval              
!               IF (r < -1d0) THEN
!                 PRINT "(A,F24.17)", "R = ", r
!                 r = -1d0
!                 CALL newton(r,dt(nd),ti,xr,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
!               ENDIF
              
!               ! Evaluate spline at specified parameter value (no distance minimixation)              
!               tpt = .5d0*dt(nd)*(r + 1d0) + ti               
!               CALL eval_cubic_spline(tpt,ti,ax(nd),bx(nd),cx(nd),dx(nd),x)
!               CALL eval_cubic_spline(tpt,ti,ay(nd),by(nd),cy(nd),dy(nd),y)              
              
              WRITE(60,*) x(j+1),y(j+1)
                       
            ENDDO  
            
            WRITE(40,"(I8,10(E24.17,1X))") n1, (x(j), j=1,ctp)
            WRITE(40,"(I8,10(E24.17,1X))") n1, (y(j), j=1,ctp)
            
            IF (i == neval-1) THEN
              WRITE(40,"(I8,10(E24.17,1X))") n2, x(ctp+1)
              WRITE(40,"(I8,10(E24.17,1X))") n2, y(ctp+1)            
            ENDIF
            
            
            PRINT*, "------------------------------------------------------------"            
            PRINT*, " "
            
          ENDDO

        ELSE
        
          WRITE(40,"(2(I8),19x,A)") 0,btype, "! Flow-specified normal flow boundary"
          WRITE(60,*) 0
              
        ENDIF 
        

      ENDDO
      CLOSE(30)
      CLOSE(40)
      CLOSE(60)
      
      
      
      
      CALL write_grid(base)

      END PROGRAM spline
      
      

      
      
