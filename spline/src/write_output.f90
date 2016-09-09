
      
      SUBROUTINE write_grid(mesh,base_name)

      USE globals, ONLY: grid,nverts
      USE grid_file_mod
      USE version, ONLY: version_information      

      IMPLICIT NONE
          
      INTEGER :: nd,el,bnd,et,nv,btype,nbseg
      TYPE(grid) :: mesh
      CHARACTER(100) :: base_name
      
      DO bnd = 1,mesh%nbou
        btype = mesh%fbseg(2,bnd)
        
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          mesh%fbseg(1,bnd) = mesh%fbseg(1,bnd) - 1
        ENDIF        
        
      ENDDO      
      
      CALL write_header("fort.14_spline",mesh%grid_name,mesh%ne,mesh%nn)  
      
      CALL write_coords(mesh%nn,mesh%xy,mesh%depth)
      
      CALL write_connectivity(mesh%ne,mesh%ect,mesh%el_type,nverts)
      
      CALL write_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds)
      
      CALL write_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds)
      
      CALL copy_footer(mesh%grid_file,"fort.14_spline")      
      
      
      OPEN(UNIT=40, FILE="fort.14_spline", POSITION="APPEND")
      CALL version_information(40)
      
      WRITE(40,"(A)") "-----------------------------------------------------------------------"           
      
      CALL write_input(40)         
      CLOSE(40)

      RETURN
      END SUBROUTINE write_grid
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE write_spline(n,bou,ndnums)
      
      USE globals, ONLY: rp,ax,bx,cx,dx,ay,by,cy,dy,dt
      
      IMPLICIT NONE
           
      INTEGER :: n,bou      
      INTEGER :: i
      INTEGER, DIMENSION(n) :: ndnums
      REAL(rp) :: t
      
      
      WRITE(30,*) n

      t = 0d0
      DO i = 1,n-1        
        WRITE(30,"(9(E25.12,1x),I7)"), ax(i,bou),bx(i,bou),cx(i,bou),dx(i,bou),ay(i,bou),by(i,bou),cy(i,bou),dy(i,bou),t,ndnums(i)
        t = t + dt(i,bou)
      ENDDO
      WRITE(30,"(9(E25.12,1x),I7)"), ax(n,bou),0d0,0d0,0d0,ay(n,bou),0d0,0d0,0d0,t,ndnums(n)

      
      
      RETURN
      END SUBROUTINE write_spline
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_cb_file(mesh)    
      
      USE globals, ONLY: grid      
      USE version, ONLY: version_information
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: bou,nd,pt
      INTEGER :: n1,n2
      INTEGER :: nbnds,btype
      INTEGER :: eind      
      CHARACTER(100) :: name
      CHARACTER(1) :: ctp_char      
      
      eind = INDEX(ADJUSTL(TRIM(mesh%grid_file)),".",.false.)   
      name = ADJUSTL(TRIM(mesh%grid_file(1:eind-1)))
      WRITE(ctp_char,"(I1)") mesh%ctp       
      
      OPEN(unit=40,file=ADJUSTL(TRIM(name)) // "_ctp" // ctp_char // ".cb")    
      
      CALL version_information(40)
      
      WRITE(40,"(A)") "-----------------------------------------------------------------------"           
      
      CALL write_input(40)           
      
      WRITE(40,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"      
      
      WRITE(40,"(I8,19x,A)") mesh%nbou, "! total number of normal flow boundaries"
      WRITE(40,"(2(I8),19x,A)") mesh%nvel,mesh%ctp, "! max number of normal flow nodes, ctp order"        
      
      DO bou = 1,mesh%nbou
      
        nbnds = mesh%fbseg(1,bou)
        btype = mesh%fbseg(2,bou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
             btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries 
             
          WRITE(40,"(2(I8),19x,A)") nbnds,btype, "! number of nodes in boundary, boundary type"               
             
          IF (nbnds > 0) THEN   
            DO nd = 1,nbnds-1
              n1 = mesh%fbnds(nd,bou)            
              WRITE(40,"(I8,1X,10(E24.17,1X))") n1, mesh%xy(1,n1), (mesh%bndxy(1,pt,nd,bou), pt=1,mesh%ctp-1)
              WRITE(40,"(I8,1X,10(E24.17,1X))") n1, mesh%xy(2,n1), (mesh%bndxy(2,pt,nd,bou), pt=1,mesh%ctp-1)             
            ENDDO
            n2 = mesh%fbnds(nbnds,bou)
            WRITE(40,"(I8,1X,10(E24.17,1X))") n2, mesh%xy(1,n2)
            WRITE(40,"(I8,1X,10(E24.17,1X))") n2, mesh%xy(2,n2) 
          
          ENDIF
        ELSE
        
          WRITE(40,"(2(I8),19x,A)") 0,btype, "! Flow-specified normal flow boundary"        
        
        ENDIF
      ENDDO
      
      CLOSE(40)
      
      RETURN
      END SUBROUTINE write_cb_file 
     


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      