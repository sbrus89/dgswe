      PROGRAM error

      USE globals
      USE read_grid
      USE evaluate
      USE basis

      IMPLICIT NONE

      INTEGER :: nd,pt,elc,elf,elb,i,dof
      INTEGER :: etf,etc,npts
      INTEGER :: ne,mndof
      INTEGER :: order
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: r,s,hb      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: xf,yf      
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hc,Qxc,Qyc
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Hf,Qxf,Qyf            
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: phi
      REAL(pres) :: HL2,QxL2,QyL2 
      REAL(pres) :: tcoarse,tfine
      REAL(pres) :: rfac

      
      CALL read_input()
      
      rfac = 2d0
      order = coarse%p + 1
      
      CALL read_grids()
      
      CALL vandermonde(coarse)
      CALL vandermonde(fine)
      CALL vandermonde(base)
            
      CALL re_vert(coarse)
      CALL re_vert(fine)
      CALL re_vert(base)
      
      CALL area_qpts()        
      CALL function_eval() ! evaluate basis and shape functions for fine grid quadrature points
      CALL detJ_eval()     ! compute determinant of Jacobian for fine grid elements at quadrature points
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      ALLOCATE(xf(mnqpta),yf(mnqpta))
      ALLOCATE(r(mnqpta),s(mnqpta),hb(mnqpta))
      ALLOCATE(elf2elc(fine%ne),elf2elb(fine%ne))
      
      CALL find_nesting(coarse,fine,elf2elc)
      CALL find_nesting(base,fine,elf2elb)

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      mndof = coarse%mndof
      ne = coarse%ne
      
      ALLOCATE(coarse%H(ne,mndof),coarse%Qx(ne,mndof),coarse%Qy(ne,mndof))
      ALLOCATE(coarse%phi(mndof,mnqpta,1))
      ALLOCATE(phi(mndof*mnqpta))
 
      mndof = fine%mndof
      ne = fine%ne
      
      ALLOCATE(fine%H(ne,mndof),fine%Qx(ne,mndof),fine%Qy(ne,mndof)) 
      
      ALLOCATE(Hc(mnqpta),Qxc(mnqpta),Qyc(mnqpta))
      ALLOCATE(Hf(mnqpta),Qxf(mnqpta),Qyf(mnqpta))
              
      CALL read_solution()
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      HL2 = 0d0
      QxL2 = 0d0
      QyL2 = 0d0
      
elemf:DO elf = 1,fine%ne        
        etf = fine%el_type(elf)
        npts = nqpta(etf)
                 
        elc = elf2elc(elf)    
        elb = elf2elb(elf)
        
        IF (fine%bndel(elf) == 1) THEN
          CYCLE elemf                
        ELSE IF(base%bndel(elb) == 1) THEN
          CYCLE elemf
        ENDIF
        
        etc = coarse%el_type(elc)  
        
        xf = 0d0
        yf = 0d0
        DO pt = 1,npts
          DO nd = 1,fine%nnds(etf)
            xf(pt) = xf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,1)
            yf(pt) = yf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,2)
          ENDDO        
        ENDDO
        
        CALL newton(coarse,xf,yf,npts,elc,r,s,hb)
                 
        IF (mod(etc,2) == 1) THEN
#ifndef adcirc        
          CALL tri_basis(coarse%p,coarse%ndof(etc),npts,r,s,phi)  
#else          
          CALL adcirc_basis(coarse%p,coarse%ndof(etc),npts,r,s,phi)    ! also need to switch in function_eval for fine solution phi
#endif
        ELSE IF (mod(etc,2) == 0) THEN
          CALL quad_basis(coarse%p,coarse%ndof(etc),npts,r,s,phi)
        ENDIF          
          
        DO pt = 1,npts
          DO dof = 1,coarse%ndof(etc)
             i = (dof-1)*npts + pt
            coarse%phi(dof,pt,1) = phi(i)
          ENDDO
        ENDDO
          
        Hc = 0d0
        Qxc = 0d0
        Qyc = 0d0
        DO pt = 1,npts
          DO dof = 1,coarse%ndof(etc)
            Hc(pt) = Hc(pt) + coarse%H(elc,dof)*coarse%phi(dof,pt,1)
            Qxc(pt) = Qxc(pt) + coarse%Qx(elc,dof)*coarse%phi(dof,pt,1)            
            Qyc(pt) = Qyc(pt) + coarse%Qy(elc,dof)*coarse%phi(dof,pt,1)
          ENDDO                
        ENDDO
        
        Hf = 0d0
        Qxf = 0d0
        Qyf = 0d0
        DO pt = 1,npts
          DO dof = 1,fine%ndof(etf)
            Hf(pt) = Hf(pt) + fine%H(elf,dof)*fine%phi(dof,pt,etf)
            Qxf(pt) = Qxf(pt) + fine%Qx(elf,dof)*fine%phi(dof,pt,etf)            
            Qyf(pt) = Qyf(pt) + fine%Qy(elf,dof)*fine%phi(dof,pt,etf)
          ENDDO                
        ENDDO       
        
        DO pt = 1,npts
          HL2 = HL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Hc(pt)-Hf(pt))**2
          QxL2 = QxL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Qxc(pt)-Qxf(pt))**2
          QyL2 = QyL2 + fine%detJ(pt,elf)*wpta(pt,etf)*(Qyc(pt)-Qyf(pt))**2
        ENDDO
        
      ENDDO elemf    
      
      HL2 = sqrt(HL2)
      QxL2 = sqrt(QxL2)
      QyL2 = sqrt(QyL2)
      
      PRINT("(A)"), "L2 Difference"
      PRINT("(A,F22.15)"), "HL2 = ", HL2
      PRINT("(A,F22.15)"), "QxL2 = ", QxL2
      PRINT("(A,F22.15)"), "QyL2 = ", QyL2
      PRINT*, " " 
      
      PRINT("(A)"), "Coarse Grid Error Estimates"
      PRINT("(A,F22.15)"), "EcH = ", HL2*rfac**order/(rfac**order-1d0)
      PRINT("(A,F22.15)"), "EcQx = ", QxL2*rfac**order/(rfac**order-1d0)
      PRINT("(A,F22.15)"), "EcQy = ", QyL2*rfac**order/(rfac**order-1d0)  
      PRINT*, " "       
      
      PRINT("(A)"), "Fine Grid Error Estimates"
      PRINT("(A,F22.15)"), "EfH = ", HL2/(rfac**order-1d0)
      PRINT("(A,F22.15)"), "EfQx = ", QxL2/(rfac**order-1d0)
      PRINT("(A,F22.15)"), "EfQy = ", QyL2/(rfac**order-1d0)         
          
      END PROGRAM error