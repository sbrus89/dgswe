      PROGRAM elongated_quad_channel

      USE globals, ONLY: rp                         
      USE grid_file_mod, ONLY: grid_type,write_grid,write_cb_file
      USE basis, ONLY: lglpts

      IMPLICIT NONE

      TYPE(grid_type) :: mesh
      INTEGER :: nd,bnd,pt,bou,seg
      INTEGER :: nbnds,btype,n,a
      INTEGER :: i,j
      INTEGER :: n1,n2
      INTEGER :: it
      REAL(rp) :: x,y,dy,dy2
      CHARACTER(100) :: forcing_file
      REAL(rp) :: Lx,Ly
      REAL(rp) :: xc,wc
      INTEGER :: nx,ny
      REAL(rp) :: ymin,ymax      
      REAL(rp) :: xe,ye
      REAL(rp) :: zmax
      REAL(rp) :: d,dx
      REAL(rp) :: th,sh
      REAL(rp) :: tol
      REAL(rp) :: flow
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r      
      CHARACTER(2) :: nx_char,ny_char
      
      Ly = 500d0
      Lx = 6000d0
      xc = 2000d0
      wc = 300d0
     
      OPEN(UNIT=10,FILE='equad.in')
      READ(10,*) nx
      READ(10,*) ny     
      READ(10,*) mesh%ctp ! ctp = 0 for adcirc
      
      WRITE(nx_char,"(I2)") nx
      WRITE(ny_char,"(I2)") ny      
      
      mesh%grid_file = "converge_equad_nx" // TRIM(ADJUSTL(nx_char)) // "_ny" // TRIM(ADJUSTL(ny_char)) // ".grd"    
      mesh%grid_name = "Elongated quad channel"
      forcing_file = "converge_equad_ny" // TRIM(ADJUSTL(ny_char)) // ".bfr"      
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create the DG grid file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      mesh%nn = nx*ny
      mesh%ne = (nx-1)*(ny-1)    
      dx = Lx/real(nx-1,rp)
      zmax = .5d0*(Ly-wc)           
      
      ALLOCATE(mesh%xy(2,mesh%nn),mesh%depth(mesh%nn))
      ALLOCATE(mesh%ect(4,mesh%ne),mesh%el_type(mesh%ne))      
      
      mesh%depth = 10d0       
      
      mesh%nope = 1
      mesh%neta = ny
      ALLOCATE(mesh%obseg(mesh%neta),mesh%obnds(ny,mesh%neta))
      mesh%obseg = 0
      
      mesh%nbou = 3
      mesh%nvel = ny + 2*nx
      ALLOCATE(mesh%fbseg(2,mesh%nbou),mesh%fbnds(mesh%nvel,mesh%nbou))
      mesh%fbseg = 0
      
      ! Calculate mesh nodes
      n = 0
      DO i = 1,nx  
        x = real(i-1,rp)*dx
        
        ymin = zmax*(1d0/(COSH(4d0*(x-xc)/Ly)))
        ymax = Ly - zmax*(1d0/(COSH(4d0*(x-xc)/Ly)))        
        dy = (ymax-ymin)/real(ny-1,rp)
        
        DO j = 1,ny          
          n = n + 1          
          
          y = ymin + real(j-1,rp)*dy 
          
          mesh%xy(1,n) = x
          mesh%xy(2,n) = y
          
          IF (i == 1) THEN
            mesh%fbseg(1,1) = mesh%fbseg(1,1) + 1
            mesh%fbseg(2,1) = 22
!             mesh%fbseg(2,1) = 12            
            mesh%fbnds(mesh%fbseg(1,1),1) = n
          ELSE IF (i == nx) THEN
            mesh%obseg(1) = mesh%obseg(1) + 1
            mesh%obnds(mesh%obseg(1),1) = n 
          ENDIF
          
          IF (j == 1) THEN
            mesh%fbseg(1,2) = mesh%fbseg(1,2) + 1
            mesh%fbseg(2,2) = 20   
!             mesh%fbseg(2,2) = 0            
            mesh%fbnds(mesh%fbseg(1,2),2) = n            
          ELSE IF (j == ny) THEN
            mesh%fbseg(1,3) = mesh%fbseg(1,3) + 1
            mesh%fbseg(2,3) = 20  
!             mesh%fbseg(2,3) = 0              
            mesh%fbnds(mesh%fbseg(1,3),3) = n             
          ENDIF
        ENDDO    
        
      ENDDO
      
      n = mesh%fbseg(1,1)
      DO i = 1,INT(n/2)
        a = mesh%fbnds(i,1)
        mesh%fbnds(i,1) = mesh%fbnds(n-i+1,1)
        mesh%fbnds(n-i+1,1) = a        
      ENDDO

      n = mesh%fbseg(1,3)
      DO i = 1,INT(n/2)
        a = mesh%fbnds(i,3)
        mesh%fbnds(i,3) = mesh%fbnds(n-i+1,3)
        mesh%fbnds(n-i+1,3) = a        
      ENDDO      
      
    
      
        ! Calculate DG element connectivity table
        n = 0     
        DO i = 1,nx-1
          DO j = 1,ny-1        
            n = n + 1
            mesh%ect(1,n) = (i-1)*ny + j
            mesh%ect(2,n) = i*ny + j
            mesh%ect(3,n) = i*ny + j + 1
            mesh%ect(4,n) = (i-1)*ny + j + 1
          ENDDO
        ENDDO
       
        mesh%el_type = 2            
      
        CALL write_grid(mesh)
       
      IF (mesh%ctp > 0) THEN         
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create the DG curved boundary file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
        ALLOCATE(mesh%bndxy(2,mesh%ctp-1,mesh%nvel,mesh%nbou))
        ALLOCATE(r(mesh%ctp+1))
      
        CALL lglpts(mesh%ctp,r)      
      
        tol = 1d-10
        DO bou = 1,mesh%nbou
      
          nbnds = mesh%fbseg(1,bou)
          btype = mesh%fbseg(2,bou)
          IF( btype == 0 .OR. btype == 10 .OR. btype == 20) THEN    ! island boundaries                          
            DO nd = 1,nbnds-1
              n1 = mesh%fbnds(nd,bou) 
              n2 = mesh%fbnds(nd+1,bou)            
              DO pt = 1,mesh%ctp-1
                xe = .5d0*((1d0-r(pt+1))*mesh%xy(1,n1) + (1d0+r(pt+1))*mesh%xy(1,n2))
                ye = .5d0*((1d0-r(pt+1))*mesh%xy(2,n1) + (1d0+r(pt+1))*mesh%xy(2,n2))              
              
                x = xe              
                DO it = 1,100
                  th = TANH(4d0*(x-xc)/Ly)
                  sh = 1d0/(COSH(4d0*(x-xc)/Ly))
              
                  IF (ye < .5d0*Ly) THEN
                    y = zmax*sh
                    dy = -zmax*th*sh*(4d0/Ly)
                    dy2 = -zmax*(1d0/(COSH(4d0*(x-xc)/Ly)**2)*sh-th**2*sh)*(4d0/Ly)**2
                  ELSE
                    y = Ly - zmax*sh     
                    dy = zmax*th*sh*(4d0*Ly)
                    dy2 = zmax*(1d0/(COSH(4d0*(x-xc)/Ly)**2)*sh-th**2*sh)*(4d0/Ly)**2
                  ENDIF                
                
                  d = 2d0*(x-xe) + 2d0*(y-ye)*dy
                  dx = 2d0 + 2d0*(dy**2 + dy2*(y-ye))
                  x = x - d/dx
                
                  IF (abs(d) < tol) THEN
                    EXIT
                  ENDIF
                ENDDO
              
                IF (ye < .5d0*Ly) THEN
                  y = zmax*(1d0/(COSH(4d0*(x-xc)/Ly)))
                ELSE
                  y = Ly - zmax*(1d0/(COSH(4d0*(x-xc)/Ly)))                   
                ENDIF
              
                mesh%bndxy(1,pt,nd,bou) = x
                mesh%bndxy(2,pt,nd,bou) = y
              ENDDO
            ENDDO
          ENDIF
        
        ENDDO     
      
        CALL write_cb_file(mesh)
        
      ENDIF
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create the DG boundary condition file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(UNIT=15, FILE=forcing_file)

      ! Write out open boundary periodic forcing data
      WRITE(15,*) 1

      WRITE(15,*) "STEADY"     
      WRITE(15,*) 0d0,1d0,0d0

      WRITE(15,*) "STEADY"
      DO seg = 1,mesh%nope
        DO nd = 1,mesh%obseg(seg)
          WRITE(15,*) 0d0,0d0
        ENDDO
      ENDDO

      ! Write out flow boundary periodic forcing data
      WRITE(15,*) 1

      WRITE(15,*) "FLOW"           
      WRITE(15,*) 0d0,1d0,0d0
      
      IF (mesh%ctp == 0) THEN
        flow = 10d0  ! double flux for adcirc criss-cross
      ELSE
        flow = 5d0
      ENDIF
      
      WRITE(15,*) "FLOW"         
      DO seg = 1,mesh%nbou
        btype = mesh%fbseg(2,seg)
        IF (btype == 2 .OR. btype == 12 .OR. btype == 22) THEN
          DO nd = 1,mesh%fbseg(1,seg)
            WRITE(15,*) flow,0d0
          ENDDO
        ENDIF
      ENDDO
     
      WRITE(15,*) 0
      
      IF (mesh%ctp == 0) THEN
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create the adcirc grid file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      DEALLOCATE(mesh%ect,mesh%el_type)
      
      mesh%ne = 4*(nx-1)*(ny-1)
      ALLOCATE(mesh%ect(3,mesh%ne),mesh%el_type(mesh%ne))
      mesh%grid_file = "converge_equad_nx" // TRIM(ADJUSTL(nx_char)) // "_ny" // TRIM(ADJUSTL(ny_char)) // "_cross.grd"         
      
      ! Calculate element connectivity table
      n = 0     
      DO i = 1,nx-1
        DO j = 1,ny-1        
          n = n + 1
          mesh%ect(1,n) = (i-1)*ny + j
          mesh%ect(2,n) = i*ny + j
          mesh%ect(3,n) = i*ny + j + 1
          
          n = n + 1
          mesh%ect(1,n) = (i-1)*ny + j
          mesh%ect(2,n) = i*ny + j + 1
          mesh%ect(3,n) = (i-1)*ny + j + 1
          
          n = n + 1
          mesh%ect(1,n) = (i-1)*ny + j
          mesh%ect(2,n) = i*ny + j
          mesh%ect(3,n) = (i-1)*ny + j + 1
          
          n = n + 1
          mesh%ect(1,n) = (i-1)*ny + j + 1
          mesh%ect(2,n) = i*ny + j
          mesh%ect(3,n) = i*ny + j + 1          
        ENDDO
      ENDDO      
      
      mesh%el_type = 1            
      
      CALL write_grid(mesh)     
      
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       ! Create the adcirc plotting grid file
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
!       DEALLOCATE(mesh%ect,mesh%el_type)
!       
!       mesh%ne = 2*(nx-1)*(ny-1)
!       ALLOCATE(mesh%ect(3,mesh%ne),mesh%el_type(mesh%ne))
!       mesh%grid_file = "converge_equad_nx" // TRIM(ADJUSTL(nx_char)) // "_ny" // TRIM(ADJUSTL(ny_char)) // "_plot.grd"           
!       
!       ! Calculate element connectivity table
!       n = 0     
!       DO i = 1,nx-1
!         DO j = 1,ny-1        
!           n = n + 1
!           mesh%ect(1,n) = (i-1)*ny + j
!           mesh%ect(2,n) = i*ny + j
!           mesh%ect(3,n) = i*ny + j + 1
!           
!           n = n + 1
!           mesh%ect(1,n) = (i-1)*ny + j
!           mesh%ect(2,n) = i*ny + j + 1
!           mesh%ect(3,n) = (i-1)*ny + j + 1
!                    
!         ENDDO
!       ENDDO      
!       
!       mesh%el_type = 1            
!       
!       CALL write_grid(mesh)   
      
      ENDIF

      END PROGRAM elongated_quad_channel
      
      
