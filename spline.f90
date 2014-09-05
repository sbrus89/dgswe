      PROGRAM spline

      USE globals, ONLY: pres,ne,nn,ctp,xy,ect,vct,depth,nelnds, &
                         nope,neta,obnds,nbou,nvel,fbseg,fbnds, &
                         nepn,epn,nverts,el_type, &
                         mnelnds

      IMPLICIT NONE
      INTEGER :: i,j,k,n,seg,sind,eind,num,qpts,btype
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert
      REAL(pres) :: htest,mult,dt,t,tpt,x,y,r,sig
      REAL(pres) :: d1,d2,d3,t1,t2
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: ax,bx,cx,dx
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: ay,by,cy,dy
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: Ml,Md,Mu,v

      OPEN(unit=30,file='spline.out')
      
      PRINT "(A)", " "

      CALL read_grid()
      
      CALL connect()
      
      sig = 1d0
      
      num = 0
      DO seg = 1,nbou
        IF(fbseg(2,seg) == 10 .OR. fbseg(2,seg) == 11)THEN
          num = num + 1
        ENDIF
      ENDDO

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A)", " "

      WRITE(30,*) num

      DO seg = 1,nbou

        IF(fbseg(2,seg) == 10 .OR. fbseg(2,seg) == 11)THEN
          n = fbseg(1,seg)    ! n nodes, n-1 subintervals

          PRINT "(A)", " "
          PRINT "(A,I5)", "Normal flow boundary ",seg
          PRINT "(A,I5)", "Normal flow boundary nodes ",n
          PRINT "(A)", " "

!           PAUSE
  

          ALLOCATE(ax(n),cx(n),bx(n-1),dx(n-1))
          ALLOCATE(ay(n),cy(n),by(n-1),dy(n-1))
          ALLOCATE(Ml(n),Md(n),Mu(n),v(n))

          dt = 1.0/(real(n,pres)-1.0)

          PRINT*, dt

          !!!!!!!!!!!!!!!!!!!
          ! x value spline 
          !!!!!!!!!!!!!!!!!!!

          ! Load nodal boundary x coordinates
          k = 1   
          DO i = 1,n
            ax(k) = xy(1,fbnds(i,seg))
            k = k+1
          ENDDO

! No tension         
!           ! Set up matrix 
!           Ml(1) = 0d0
!           Md(1) = 1d0
!           Mu(1) = 0d0
!           DO i = 2,n-1
!             Ml(i) = dt
!             Md(i) = 4d0*dt
!             Mu(i) = dt
!           ENDDO
!           Ml(n) = 0d0
!           Md(n) = 1d0
!           Mu(n) = 0d0
! 
!           ! Set up RHS
!           v(1) = 0d0
!           DO i = 2,n-1
!             v(i) = 3d0/dt*(ax(i+1)-2d0*ax(i)+ax(i-1))
!           ENDDO
!           v(n) = 0d0

! With tension
          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1
            Ml(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
            Md(i) = (4d0/sig**2)*((sig*cosh(sig*dt))/sinh(sig*dt)-1d0/dt)
            Mu(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          v(1) = 0d0
          DO i = 2,n-1
            v(i) = 1d0/dt*(ax(i+1)-2d0*ax(i)+ax(i-1))
          ENDDO
          v(n) = 0d0

          ! Solve system for c coefficients, forward sweep
          DO i = 2,n
            mult = Ml(i)/Md(i-1)
            Md(i) = Md(i) - mult*Mu(i-1)
            v(i) = v(i) - mult*v(i-1)
          ENDDO

          ! Solve system for c coefficients, backward sweep
          cx(n) = v(n)/Md(n)
          DO i = n-1,1,-1
            cx(i) = (v(i) - Mu(i)*cx(i+1))/Md(i)
          ENDDO

          ! Solve for other coefficients d and b
          DO i = 1,n-1
            dx(i) = (cx(i+1)-cx(i))/(3.0*dt)
            bx(i) = (ax(i+1)-ax(i))/dt - dt*(2.0*cx(i)+cx(i+1))/3.0
          ENDDO

          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!

          ! Load nodal boundary y coordinates
          k = 1
          DO i = 1,n
            ay(k) = xy(2,fbnds(i,seg))
            k = k+1
          ENDDO

! No tension          
!           ! Set up matrix 
!           Ml(1) = 0d0
!           Md(1) = 1d0
!           Mu(1) = 0d0
!           DO i = 2,n-1
!             Ml(i) = dt
!             Md(i) = 4d0*dt
!             Mu(i) = dt
!           ENDDO
!           Ml(n) = 0d0
!           Md(n) = 1d0
!           Mu(n) = 0d0
! 
!           ! Set up RHS
!           v(1) = 0d0
!           DO i = 2,n-1
!             v(i) = 3d0/dt*(ay(i+1)-2d0*ay(i)+ay(i-1))
!           ENDDO
!           v(n) = 0d0
          
! With tension          
          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1
            Ml(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
            Md(i) = (4d0/sig**2)*((sig*cosh(sig*dt))/sinh(sig*dt)-1d0/dt)
            Mu(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0       
          
          ! Set up RHS
          v(1) = 0d0
          DO i = 2,n-1
            v(i) = 1d0/dt*(ay(i+1)-2d0*ay(i)+ay(i-1))
          ENDDO
          v(n) = 0d0          

          ! Solve system for c coefficients, forward sweep
          DO i = 2,n
            mult = Ml(i)/Md(i-1)
            Md(i) = Md(i) - mult*Mu(i-1)
            v(i) = v(i) - mult*v(i-1)
          ENDDO

          ! Solve system for c coefficients, backward sweep
          cy(n) = v(n)/Md(n)
          DO i = n-1,1,-1
            cy(i) = (v(i) - Mu(i)*cy(i+1))/Md(i)
          ENDDO

          ! Solve for other coefficients d and b
          DO i = 1,n-1
            dy(i) = (cy(i+1)-cy(i))/(3.0*dt)
            by(i) = (ay(i+1)-ay(i))/dt - dt*(2.0*cy(i)+cy(i+1))/3.0
          ENDDO
          
          
          t = 0d0
          DO nd = 1,n-1
            n1bed = fbnds(nd,seg)
            n2bed = fbnds(nd+1,seg)    

      elem: DO el = 1,nepn(n1bed)
              eln = epn(el,n1bed)
              
              nvert = nverts(el_type(eln))
              
              DO led = 1,nvert
 
                n1ed1 = vct(mod(led+0,nvert)+1,eln)
                n2ed1 = vct(mod(led+1,nvert)+1,eln) 

                IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
                   ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
                   
                   print*, "  ", eln
                   
                   DO i = 1,ctp-1
                     r = -1d0 + real(i,pres)*2d0/real(ctp,pres)
                     tpt = .5d0*dt*(r + 1d0) + t
                     
                     x = ax(nd) + bx(nd)*(tpt-t) + cx(nd)*(tpt-t)**2 + dx(nd)*(tpt-t)**3
                     y = ay(nd) + by(nd)*(tpt-t) + cy(nd)*(tpt-t)**2 + dy(nd)*(tpt-t)**3
                     
                     ndn = mod(led,nvert)*ctp+i+1
                     
                     xy(1,ect(ndn,eln)) = x
                     xy(2,ect(ndn,eln)) = y
                   ENDDO
                   
                   t = t + dt
                   EXIT elem

                ENDIF
                
              ENDDO
            ENDDO elem
          ENDDO          

  
          WRITE(30,*) n,dt

          DO i = 1,n-1
            WRITE(30,"(4(E25.12))"), ax(i),bx(i),cx(i),dx(i)
          ENDDO
          WRITE(30,"(4(E25.12))"), ax(n),0.0,0.0,0.0

          DO i = 1,n-1
            WRITE(30,"(4(E25.12))"), ay(i),by(i),cy(i),dy(i)
          ENDDO
          WRITE(30,"(4(E25.12))"), ay(n),0.0,0.0,0.0


          DEALLOCATE(ax,bx,cx,dx)
          DEALLOCATE(ay,by,cy,dy)
          DEALLOCATE(Ml,Md,Mu,v)

        ENDIF

      ENDDO
      CLOSE(30)
      
      CALL write_grid()

      END PROGRAM spline