      MODULE transformation
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      

      CONTAINS
      
      SUBROUTINE init_element_coordinates(ne,mnnds,el_type,nverts,xy,ect,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: mnnds
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: elxy
            
      INTEGER :: el,nd
      INTEGER :: et,nv
      INTEGER :: alloc_status
      
      ALLOCATE(elxy(mnnds,ne,2), STAT=alloc_status)
      
      
      DO el = 1,ne        
        et = el_type(el)
        nv = nverts(et)
        DO nd = 1,nv
          elxy(nd,el,1) = xy(1,ect(nd,el))
          elxy(nd,el,2) = xy(2,ect(nd,el))
        ENDDO              
      ENDDO          
      
      RETURN
      END SUBROUTINE init_element_coordinates
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE element_transformation(nnd,x,y,l,xpt,ypt,dldr,dlds,drdx,drdy,dsdx,dsdy,detJ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nnd

      REAL(rp), DIMENSION(:), INTENT(IN) :: x,y
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: xpt,ypt
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(OUT), OPTIONAL :: drdx,drdy,dsdx,dsdy,detJ

      INTEGER :: nd
      REAL(rp) :: dxdr,dxds,dydr,dyds
      INTEGER :: calc_inv


      IF ( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
           PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
           PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. &
           PRESENT(detJ) ) THEN
        calc_inv = 1
      ELSE
        calc_inv = 0
      ENDIF
     

      xpt = 0d0
      ypt = 0d0
          
      DO nd = 1,nnd        
              
        xpt = xpt + l(nd)*x(nd)
        ypt = ypt + l(nd)*y(nd)
                      
      ENDDO

         
      IF (calc_inv == 1) THEN

        dxdr = 0d0
        dxds = 0d0
        dydr = 0d0
        dyds = 0d0
            
          
        DO nd = 1,nnd
        
          dxdr = dxdr + dldr(nd)*x(nd)
          dxds = dxds + dlds(nd)*x(nd)
          dydr = dydr + dldr(nd)*y(nd)
          dyds = dyds + dlds(nd)*y(nd)              
                      
        ENDDO
            
        detJ = dxdr*dyds - dxds*dydr
            
        drdx =  dyds/detJ
        drdy = -dxds/detJ
        dsdx = -dydr/detJ
        dsdy =  dxdr/detJ

      ENDIF

      END SUBROUTINE element_transformation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE xy2rs(nv,p,elx,ely,npt,x,y,r,s)

      USE shape_functions_mod, ONLY: shape_functions_area_eval

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nv
      INTEGER, INTENT(IN) :: p
      REAL(rp), DIMENSION(:), INTENT(IN) :: elx
      REAL(rp), DIMENSION(:), INTENT(IN) :: ely
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: x   
      REAL(rp), DIMENSION(:), INTENT(IN) :: y       
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r   
      REAL(rp), DIMENSION(:), INTENT(OUT) :: s      
      
      INTEGER :: it,n,pt
      INTEGER :: maxit
      INTEGER :: mnnds
      REAL(rp) :: tol    
      REAL(rp) :: f,g
      REAL(rp) :: error(npt),err
      REAL(rp) :: drdf,drdg,dsdf,dsdg,jac
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: l,dldr,dlds
        
      tol = 1d-9
      maxit = 100
          
      mnnds = (p+1)**2      
      ALLOCATE(l(mnnds,npt),dldr(mnnds,npt),dlds(mnnds,npt))
        
        
      IF (mod(nv,2) == 1) THEN
        DO pt = 1,npt
          r(pt) = -1d0/3d0
          s(pt) = -1d0/3d0
        ENDDO
      ELSE IF (mod(nv,2) == 0) THEN
        DO pt = 1,npt
          r(pt) = 0d0
          s(pt) = 0d0
        ENDDO
      ENDIF
      

      DO it = 1,maxit     
           

        CALL shape_functions_area_eval(nv,p,n,npt,r,s,l,dldr,dlds)     
        
        DO pt = 1,npt
          CALL element_transformation(n,elx,ely,l(:,pt),f,g, &
                                    dldr(:,pt),dlds(:,pt),drdf,drdg,dsdf,dsdg,jac)               
        
          f = f - x(pt)
          g = g - y(pt)
        
          r(pt) = r(pt) - (drdf*f + drdg*g)
          s(pt) = s(pt) - (dsdf*f + dsdg*g)
          
          error(pt) = MAX(ABS(f),ABS(g))
        ENDDO
              
        err = MAXVAL(error)
        IF ( err < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
!       IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",err
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",err
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ENDIF


      RETURN
      END SUBROUTINE xy2rs
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      END MODULE transformation