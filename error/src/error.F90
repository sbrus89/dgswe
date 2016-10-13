      PROGRAM error

      USE globals
      USE allocation
      USE evaluate
      USE basis
      USE area_qpts_mod
      USE transformation, ONLY: init_vandermonde,xy2rs
      USE version, ONLY: version_information

      IMPLICIT NONE

      INTEGER :: el,nd,pt,elc,elf,elb,i,dof
      INTEGER :: et,etf,etc,npts
      INTEGER :: ne,mndof
      INTEGER :: myrank,p_int,ctp_int 

      REAL(rp) :: HL2,QxL2,QyL2 
      REAL(rp) :: tcoarse,tfine
      REAL(rp) :: rfac,order
      REAL(rp) :: Hdiff,Qxdiff,Qydiff
      REAL(rp) :: Hmax,Qxmax,Qymax
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Hel,Qxel,Qyel
      
      myrank = 0
      p_int = 5
      ctp_int = 0

      CALL version_information(unit=6) ! print out current git branch/SHA
      
      CALL read_error_input()  ! read error.inp file
      
      rfac = 2d0
      order = real(coarse%p,rp) + 1d0
      
      CALL read_grids()  ! read coarse, fine, and base grids
      
      CALL connect(base)
      CALL connect(fine)      
      CALL connect(coarse)                   
      
      CALL area_qpts(myrank,p_int,ctp_int,nel_type,nqpta,mnqpta,wpta,qpta)    ! get area quadrature points   
      CALL function_eval()                                                    ! evaluate basis and shape functions for fine grid quadrature points
      CALL detJ_eval()                                                        ! compute determinant of Jacobian for fine grid elements at quadrature points
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      CALL nest_alloc()      
      CALL find_nesting(coarse,fine,elf2elc)   ! determine element nesting between coarse and fine grids
      CALL find_nesting(base,fine,elf2elb)     ! determine element nesting between base and fine grids

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL error_alloc()              
      CALL read_solution(coarse,tcoarse)  ! read in solutions    
      CALL read_solution(fine,tfine)
      
      IF(tcoarse /= tfine) THEN
        PRINT("(A)"), "Warning: time snaps do not agree"
        PRINT("(A,F15.7)"), "tcoarse = ", tcoarse
        PRINT("(A,F15.7)"), "tfine = ", tfine
        PRINT*, " " 
      ENDIF      
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "        Calculating Error Integral           "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      CALL init_vandermonde(nel_type,coarse%np)
      
      HL2 = 0d0
      QxL2 = 0d0
      QyL2 = 0d0
      
elemf:DO elf = 1,fine%ne  ! Calculate error integral in the fine grid elements

        etf = fine%el_type(elf) ! element type for fine element
        npts = nqpta(etf)       ! number of quadrature points needed
                 
        elc = elf2elc(elf) ! find element from coarse grid
        elb = elf2elb(elf) ! find element from base grid
        
        IF (fine%bndel(elf) == 1) THEN     ! Ignore elements on fine boundary
          CYCLE elemf                
        ELSE IF(base%bndel(elb) == 1) THEN ! Ignore element inside a base grid boundary element
          CYCLE elemf
        ENDIF
        
        etc = coarse%el_type(elc)  ! element type for coarse element
        
        ! evaluate x,y coordinates of fine element quadrature points
        xf = 0d0
        yf = 0d0
        DO pt = 1,npts
          DO nd = 1,fine%nnds(etf)
            xf(pt) = xf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,1)  
            yf(pt) = yf(pt) + fine%l(nd,pt,etf)*fine%elxy(nd,elf,2)
          ENDDO        
        ENDDO
        
        ! calculate r,s coordinates of fine element quadrature points
!         CALL newton(coarse,xf,yf,npts,elc,r,s,hb)  
        CALL xy2rs(etc,coarse%np(etc),coarse%elxy(:,elc,1),coarse%elxy(:,elc,2),npts,xf,yf,r,s)
               
        ! evalute coarse element basis functions at fine element quadrautre points               
        IF (mod(etc,2) == 1) THEN 
#ifdef dgswem  
          CALL dgswem_basis(coarse%p,coarse%ndof(etc),npts,r,s,coarse%phi(:,:,1))  
#elif adcirc         
          CALL linear_basis(npts,r,s,coarse%phi(:,:,1))     
#else          
          CALL tri_basis(coarse%p,coarse%ndof(etc),npts,r,s,coarse%phi(:,:,1))            
#endif
        ELSE IF (mod(etc,2) == 0) THEN
          CALL quad_basis(coarse%p,coarse%ndof(etc),npts,r,s,coarse%phi(:,:,1))
        ENDIF          
          

          
        ! evaluate coarse element solution at quadrature points  
        Hc = 0d0
        Qxc = 0d0
        Qyc = 0d0
        DO pt = 1,npts
          DO dof = 1,coarse%ndof(etc)
            Hc(pt) = Hc(pt) + coarse%H(elc,dof,1)*coarse%phi(dof,pt,1)
            Qxc(pt) = Qxc(pt) + coarse%Qx(elc,dof,1)*coarse%phi(dof,pt,1)            
            Qyc(pt) = Qyc(pt) + coarse%Qy(elc,dof,1)*coarse%phi(dof,pt,1)
          ENDDO                
        ENDDO
        
        ! evaluate fine element solution at quadrature points
        Hf = 0d0
        Qxf = 0d0
        Qyf = 0d0
        DO pt = 1,npts
          DO dof = 1,fine%ndof(etf)
            Hf(pt) = Hf(pt) + fine%H(elf,dof,1)*fine%phi(dof,pt,etf)
            Qxf(pt) = Qxf(pt) + fine%Qx(elf,dof,1)*fine%phi(dof,pt,etf)            
            Qyf(pt) = Qyf(pt) + fine%Qy(elf,dof,1)*fine%phi(dof,pt,etf)
          ENDDO                
        ENDDO       
        
        ! calcuate L2 error integral
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
      PRINT("(A,F22.15)"), "HL2  = ", HL2      
      PRINT("(A,ES22.15)"), "     = ", HL2       
      PRINT("(A,F22.15)"), "QxL2 = ", QxL2
      PRINT("(A,ES22.15)"), "     = ", QxL2      
      PRINT("(A,F22.15)"), "QyL2 = ", QyL2
      PRINT("(A,ES22.15)"), "     = ", QyL2      
      PRINT*, " " 
      
      PRINT("(A)"), "Coarse Grid Error Estimates"
      PRINT("(A,F22.15)"), "EcH  = ", HL2*rfac**order/(rfac**order-1d0)
      PRINT("(A,ES22.15)"), "     = ", HL2*rfac**order/(rfac**order-1d0)      
      PRINT("(A,F22.15)"), "EcQx = ", QxL2*rfac**order/(rfac**order-1d0)
      PRINT("(A,ES22.15)"), "     = ", QxL2*rfac**order/(rfac**order-1d0)      
      PRINT("(A,F22.15)"), "EcQy = ", QyL2*rfac**order/(rfac**order-1d0)  
      PRINT("(A,ES22.15)"), "     = ", QyL2*rfac**order/(rfac**order-1d0)        
      PRINT*, " "       
      
      PRINT("(A)"), "Fine Grid Error Estimates"
      PRINT("(A,F22.15)"), "EfH  = ", HL2/(rfac**order-1d0)
      PRINT("(A,ES22.15)"), "     = ", HL2/(rfac**order-1d0)      
      PRINT("(A,F22.15)"), "EfQx = ", QxL2/(rfac**order-1d0)
      PRINT("(A,ES22.15)"), "     = ", QxL2/(rfac**order-1d0)      
      PRINT("(A,F22.15)"), "EfQy = ", QyL2/(rfac**order-1d0)
      PRINT("(A,ES22.15)"), "     = ", QyL2/(rfac**order-1d0)       
      PRINT*, " "             
      

      
      IF(coarse%grid_file == fine%grid_file .and. coarse%p == fine%p) THEN
      
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "       Calculating Max DOF Difference        "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "            
      
        Hmax = 0d0
        Qxmax = 0d0
        Qymax = 0d0
        
        ALLOCATE(Hel(fine%ne),Qxel(fine%ne),Qyel(fine%ne))

        Hel(:) = 0d0
        Qxel(:) = 0d0
        Qyel(:) = 0d0
        
        DO el = 1,fine%ne
          et = fine%el_type(el)               
        
          DO dof = 1,coarse%ndof(et)
          
            Hdiff =  abs(coarse%H(el,dof,1)-fine%H(el,dof,1))  
            IF (Hdiff > Hmax) THEN 
              Hmax = Hdiff
            ENDIF 
            IF (Hdiff > Hel(el)) THEN
              Hel(el) = Hdiff
            ENDIF
            
            Qxdiff = abs(coarse%Qx(el,dof,1)-fine%Qx(el,dof,1))          
            IF (Qxdiff > Qxmax) THEN
              Qxmax = Qxdiff
            ENDIF
            IF (Qxdiff > Qxel(el)) THEN
              Qxel(el) = Qxdiff
            ENDIF
            
            Qydiff = abs(coarse%Qy(el,dof,1)-fine%Qy(el,dof,1))  
            IF (Qydiff > Qymax) THEN
              Qymax = Qydiff
            ENDIF
            IF (Qydiff > Qyel(el)) THEN
              Qyel(el) = Qydiff
            ENDIF
            
          ENDDO          
        ENDDO
        
        PRINT("(A)"), "Max DOF Differences"
        PRINT("(A,ES22.15)"), "Hmax = ", Hmax
        PRINT("(A,ES22.15)"), "Qxmax = ", Qxmax
        PRINT("(A,ES22.15)"), "Qymax = ", Qymax
        PRINT*, " " 
        
        OPEN(UNIT=10,FILE='maxdiff_el.d')
        
        WRITE(10,"(I8)") fine%ne
        DO el = 1,fine%ne
          WRITE(10,"(3(E24.17))") Hel(el),Qxel(el),Qyel(el)
        ENDDO
        
        CLOSE(10)
        
      ENDIF
          
      END PROGRAM error