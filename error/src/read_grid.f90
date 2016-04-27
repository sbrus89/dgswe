
      SUBROUTINE read_grids()
      
      USE globals, ONLY: coarse,fine,base      
      USE allocation, ONLY: sizes
      
      IMPLICIT NONE
      
      INTEGER :: i,j
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Coarse Grid Information           "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      CALL sizes(coarse)
      CALL read_grid(coarse)      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Fine Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
     
      CALL sizes(fine)
      CALL read_grid(fine)
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Base Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
      
      CALL sizes(base)
      CALL read_grid(base)      
      
      RETURN
      END SUBROUTINE read_grids
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,solution,nverts
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation

      IMPLICIT NONE

                 
      INTEGER :: el,j           
      INTEGER :: et,nv                        
      INTEGER :: coord_sys
      REAL(rp) :: h0
      INTEGER :: alloc_status   
      TYPE(solution) :: mesh      
      
      coord_sys = 1
      h0 = 0d0
      
      CALL read_header(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)
      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
!       CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)      
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
      
      
      
      
      
      
      ALLOCATE(mesh%elxy(mesh%mnnds,mesh%ne,2))
      ALLOCATE(mesh%elhb(mesh%mnnds,mesh%ne))

      DO el = 1,mesh%ne       
        et = mesh%el_type(el)
        nv = nverts(et)
        DO j = 1,nv
          mesh%elxy(j,el,1) = mesh%xy(1,mesh%ect(j,el))
          mesh%elxy(j,el,2) = mesh%xy(2,mesh%ect(j,el))
          mesh%elhb(j,el)   = mesh%depth(mesh%ect(j,el))
        ENDDO      
        
      ENDDO
      
      

      

      

      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds) 
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)      
      
         

      
!       ! Find elements on land boundaries (these are left out of the L2 error calculation in case they're curve)
!       ALLOCATE(sol%bndel(sol%ne))
!       sol%bndel = 0
!       
!       IF (exclude_bndel) THEN
!       DO i = 1,sol%nbou
!         btype = sol%fbseg(2,i)
!         nbseg = sol%fbseg(1,i)-1
!         IF (btype == 0 .OR. btype == 10 .OR. btype == 20 .OR. &
!             btype == 2 .OR. btype == 12 .OR. btype == 22) THEN
!           DO j = 1,nbseg
!             n1 = sol%fbnds(j,i)
!             n2 = sol%fbnds(j+1,i)
!       elem: DO el = 1,sol%ne
!               found = 0
!               DO nd = 1,sol%nelnds(el)
!                 IF (sol%ect(nd,el) == n1 .OR. sol%ect(nd,el) == n2) THEN
!                   found = found + 1
!                 ENDIF
!               ENDDO
!               IF (found == 2) THEN
!                 sol%bndel(el) = 1
!                 EXIT elem
!               ENDIF
!             ENDDO elem
!           ENDDO
!         ENDIF
!       ENDDO
! 
!       DO i = 1,sol%nope
!         nbseg = sol%obseg(i)-1
!           DO j = 1,nbseg
!             n1 = sol%obnds(j,i)
!             n2 = sol%obnds(j+1,i)
!      elem2: DO el = 1,sol%ne
!               found = 0
!               DO nd = 1,sol%nelnds(el)
!                 IF (sol%ect(nd,el) == n1 .OR. sol%ect(nd,el) == n2) THEN
!                   found = found + 1
!                 ENDIF
!               ENDDO
!               IF (found == 2) THEN
!                 sol%bndel(el) = 1
!                 EXIT elem2
!               ENDIF
!             ENDDO elem2
!           ENDDO
!       ENDDO
!       ENDIF


      RETURN
      END SUBROUTINE read_grid