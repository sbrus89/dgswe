      SUBROUTINE area_qpts( )

      USE globals, ONLY: pres,nel_type,wpta,qpta,nqpta,mnqpta

      IMPLICIT NONE 
      INTEGER :: i,pt,et
      INTEGER :: order(nel_type)
      INTEGER :: npt(nel_type)
      REAL(pres) :: w(13**2,nel_type),r(13**2,2,nel_type)
      
      order(1) = 10
      order(2) = order(1)
      order(3) = 10
      order(4) = order(3)
      
      et = 1
      DO i = 1,nel_type
      
        CALL cubature(et,order(i),npt(i),w(:,i),r(:,:,i))
        et = et*-1
        
      ENDDO
      
      
      mnqpta = maxval(npt)
      
      
      ALLOCATE(wpta(mnqpta,nel_type),qpta(mnqpta,2,nel_type))
      
      DO i = 1,nel_type
        nqpta(i) = npt(i)
        DO pt = 1,mnqpta
          wpta(pt,i) = w(pt,i)
          qpta(pt,1,i) = r(pt,1,i)
          qpta(pt,2,i) = r(pt,2,i)
        ENDDO
      ENDDO          


      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "        Area Integration Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

!       PRINT "(A)", 'Area quadrature weights and points (r,s)'
      
      DO i = 1,4
      PRINT "(A,I5)", "Number of area quadrature points:",nqpta(i)      
!         DO pt = 1,nqpta(i)
!           PRINT "(3(F10.3))", wpta(pt,i),qpta(pt,1,i),qpta(pt,2,i)
!         ENDDO
!         PRINT*, ' '
      ENDDO
      PRINT*, " "

      
      RETURN
      END SUBROUTINE area_qpts
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE cubature(eltype,p,nqpta,wpta,qpta)
      
      USE globals, ONLY: pres
      
      IMPLICIT NONE
      
      INTEGER :: p,eltype
      INTEGER :: i,j,n
      INTEGER :: nqpta,npt
      REAL(pres) :: wpta(13**2),qpta(13**2,2)
      REAL(pres) :: w(13),r(13)
      
      
      wpta = 0d0
      qpta = 0d0
      
      IF (eltype == 1) THEN
      
        SELECT CASE( p ) 
          CASE( 1 )
            nqpta  = 1 
          
            qpta( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01  /)  
            wpta(1) = 2.000000000000000d+00

          CASE( 2 )
            nqpta = 3

            qpta( 1, :) = (/ -6.666666666666670d-01,  -6.666666666666670d-01  /) 
            qpta( 2, :) = (/ 3.333333333333330d-01,  -6.666666666666670d-01  /) 
            qpta( 3, :) = (/ -6.666666666666670d-01,  3.333333333333330d-01  /)  
            wpta(1) = 6.666666666666670d-01
            wpta(2) = 6.666666666666670d-01
            wpta(3) = 6.666666666666670d-01

          CASE( 3 )
            nqpta = 6

            qpta( 1, :) = (/ -8.168475729804580d-01,  -8.168475729804580d-01  /)  
            qpta( 2, :) = (/ 6.336951459609170d-01,  -8.168475729804590d-01  /)  
            qpta( 3, :) = (/ -8.168475729804590d-01,  6.336951459609170d-01  /)  
            qpta( 4, :) = (/ -1.081030181680700d-01,  -1.081030181680700d-01  /)  
            qpta( 5, :) = (/ -7.837939636638600d-01,  -1.081030181680700d-01  /)  
            qpta( 6, :) = (/ -1.081030181680700d-01,  -7.837939636638600d-01  /)  
            wpta(1) = 2.199034873106440d-01
            wpta(2) = 2.199034873106440d-01
            wpta(3) = 2.199034873106440d-01
            wpta(4) = 4.467631793560230d-01
            wpta(5) = 4.467631793560230d-01
            wpta(6) = 4.467631793560230d-01

          CASE( 4 )
            nqpta = 6

            qpta( 1, :) = (/ -8.168475729804580d-01,  -8.168475729804580d-01  /)  
            qpta( 2, :) = (/ 6.336951459609170d-01,  -8.168475729804590d-01  /)  
            qpta( 3, :) = (/ -8.168475729804590d-01,  6.336951459609170d-01  /)  
            qpta( 4, :) = (/ -1.081030181680700d-01,  -1.081030181680700d-01  /)  
            qpta( 5, :) = (/ -7.837939636638600d-01,  -1.081030181680700d-01  /)  
            qpta( 6, :) = (/ -1.081030181680700d-01,  -7.837939636638600d-01  /)  
            wpta(1) = 2.199034873106440d-01
            wpta(2) = 2.199034873106440d-01
            wpta(3) = 2.199034873106440d-01
            wpta(4) = 4.467631793560230d-01
            wpta(5) = 4.467631793560230d-01
            wpta(6) = 4.467631793560230d-01 

          CASE( 5 )
            nqpta = 7
  
            qpta( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01  /)  
            qpta( 2, :) = (/ -5.971587178977000d-02,  -5.971587178977000d-02  /)  
            qpta( 3, :) = (/ -8.805682564204600d-01,  -5.971587178977000d-02  /)  
            qpta( 4, :) = (/ -5.971587178977000d-02,  -8.805682564204600d-01  /)  
            qpta( 5, :) = (/ -7.974269853530870d-01,  -7.974269853530870d-01  /)  
            qpta( 6, :) = (/ 5.948539707061750d-01,  -7.974269853530870d-01  /)  
            qpta( 7, :) = (/ -7.974269853530870d-01,  5.948539707061750d-01  /)  
            wpta(1) = 4.500000000000000d-01
            wpta(2) = 2.647883055770120d-01
            wpta(3) = 2.647883055770120d-01
            wpta(4) = 2.647883055770120d-01
            wpta(5) = 2.518783610896540d-01
            wpta(6) = 2.518783610896540d-01
            wpta(7) = 2.518783610896540d-01

          CASE( 6 )
            nqpta = 12

            qpta( 1, :) = (/ -5.014265096581790d-01,  -5.014265096581790d-01  /)  
            qpta( 2, :) = (/ 2.853019316358000d-03,  -5.014265096581790d-01  /)  
            qpta( 3, :) = (/ -5.014265096581790d-01,  2.853019316358000d-03  /)  
            qpta( 4, :) = (/ -8.738219710169960d-01,  -8.738219710169960d-01  /)  
            qpta( 5, :) = (/ 7.476439420339910d-01,  -8.738219710169960d-01  /)  
            qpta( 6, :) = (/ -8.738219710169960d-01,  7.476439420339910d-01  /)  
            qpta( 7, :) = (/ -3.792950979324310d-01,  -8.937099003103660d-01  /)  
            qpta( 8, :) = (/ -8.937099003103660d-01,  -3.792950979324310d-01  /)  
            qpta( 9, :) = (/ 2.730049982427970d-01,  -8.937099003103660d-01  /)  
            qpta( 10, :) = (/ -8.937099003103660d-01,  2.730049982427970d-01  /)  
            qpta( 11, :) = (/ 2.730049982427970d-01,  -3.792950979324310d-01  /)  
            qpta( 12, :) = (/ -3.792950979324310d-01,  2.730049982427970d-01  /) 
            wpta(1) = 2.335725514527590d-01
            wpta(2) = 2.335725514527590d-01
            wpta(3) = 2.335725514527590d-01
            wpta(4) = 1.016898127404140d-01
            wpta(5) = 1.016898127404140d-01
            wpta(6) = 1.016898127404140d-01
            wpta(7) = 1.657021512367470d-01 
            wpta(8) = 1.657021512367470d-01
            wpta(9) = 1.657021512367470d-01
            wpta(10) = 1.657021512367470d-01
            wpta(11) = 1.657021512367470d-01
            wpta(12) = 1.657021512367470d-01

          CASE( 7 )
            nqpta = 15

            qpta( 1, :) = (/ -1.587650707239070d-01,  -1.587650707239070d-01  /)  
            qpta( 2, :) = (/ -6.824698585521860d-01,  -1.587650707239070d-01  /)  
            qpta( 3, :) = (/ -1.587650707239070d-01,  -6.824698585521860d-01  /)  
            qpta( 4, :) = (/ -9.019376324216080d-01,  -9.019376324216080d-01  /)  
            qpta( 5, :) = (/ 8.038752648432150d-01,  -9.019376324216080d-01  /)  
            qpta( 6, :) = (/ -9.019376324216080d-01,  8.038752648432150d-01  /)  
            qpta( 7, :) = (/ -6.966960025905280d-01,  -6.966960025905280d-01  /)  
            qpta( 8, :) = (/ 3.933920051810560d-01,  -6.966960025905280d-01  /)  
            qpta( 9, :) = (/ -6.966960025905280d-01,  3.933920051810560d-01  /)  
            qpta( 10, :) = (/ 3.482755759622830d-01,  -3.757620558562780d-01  /)  
            qpta( 11, :) = (/ -3.757620558562780d-01,  3.482755759622830d-01  /)  
            qpta( 12, :) = (/ -9.725135201060050d-01,  -3.757620558562780d-01  /)  
            qpta( 13, :) = (/ -3.757620558562780d-01,  -9.725135201060050d-01  /)  
            qpta( 14, :) = (/ -9.725135201060050d-01,  3.482755759622830d-01  /)  
            qpta( 15, :) = (/ 3.482755759622830d-01,  -9.725135201060050d-01  /)  
            wpta(1) = 2.784029593504940d-01
            wpta(2) = 2.784029593504940d-01
            wpta(3) = 2.784029593504940d-01
            wpta(4) = 6.365388265714500d-02
            wpta(5) = 6.365388265714500d-02
            wpta(6) = 6.365388265714500d-02
            wpta(7) = 1.686426467597030d-01 
            wpta(8) = 1.686426467597030d-01
            wpta(9) = 1.686426467597030d-01
            wpta(10) = 7.798358894966199d-02
            wpta(11) = 7.798358894966199d-02
            wpta(12) = 7.798358894966199d-02
            wpta(13) = 7.798358894966199d-02
            wpta(14) = 7.798358894966199d-02
            wpta(15) = 7.798358894966199d-02

          CASE( 8 )
            nqpta = 16

            qpta( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01  /)  
            qpta( 2, :) = (/ -8.141482341455400d-02,  -8.141482341455400d-02  /)  
            qpta( 3, :) = (/ -8.371703531708929d-01,  -8.141482341455400d-02  /)  
            qpta( 4, :) = (/ -8.141482341455400d-02,  -8.371703531708929d-01  /)  
            qpta( 5, :) = (/ -6.588613844964800d-01,  -6.588613844964800d-01  /)  
            qpta( 6, :) = (/ 3.177227689929590d-01,  -6.588613844964800d-01  /)  
            qpta( 7, :) = (/ -6.588613844964800d-01,  3.177227689929590d-01  /)  
            qpta( 8, :) = (/ -8.989055433659380d-01,  -8.989055433659380d-01  /)  
            qpta( 9, :) = (/ 7.978110867318760d-01,  -8.989055433659380d-01  /)  
            qpta( 10, :) = (/ -8.989055433659380d-01,  7.978110867318760d-01  /)  
            qpta( 11, :) = (/ -4.737743407307240d-01,  4.569847859108090d-01  /)  
            qpta( 12, :) = (/ 4.569847859108090d-01,  -4.737743407307240d-01  /)  
            qpta( 13, :) = (/ -9.832104451800850d-01,  4.569847859108090d-01  /)  
            qpta( 14, :) = (/ 4.569847859108090d-01,  -9.832104451800850d-01  /)  
            qpta( 15, :) = (/ -9.832104451800850d-01,  -4.737743407307240d-01  /)  
            qpta( 16, :) = (/ -4.737743407307240d-01,  -9.832104451800850d-01  /)
            wpta(1) = 2.886312153555740d-01
            wpta(2) = 1.901832685345690d-01
            wpta(3) = 1.901832685345690d-01
            wpta(4) = 1.901832685345690d-01
            wpta(5) = 2.064347410694370d-01
            wpta(6) = 2.064347410694370d-01
            wpta(7) = 2.064347410694370d-01 
            wpta(8) = 6.491699524639601d-02
            wpta(9) = 6.491699524639601d-02
            wpta(10) = 6.491699524639601d-02
            wpta(11) = 5.446062834887000d-02
            wpta(12) = 5.446062834887000d-02
            wpta(13) = 5.446062834887000d-02
            wpta(14) = 5.446062834887000d-02
            wpta(15) = 5.446062834887000d-02
            wpta(16) = 5.446062834887000d-02
  
          CASE( 9 )
            nqpta = 19
          
            qpta( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01  /)  
            qpta( 2, :) = (/ -2.063496160252500d-02,  -2.063496160252500d-02  /)  
            qpta( 3, :) = (/ -9.587300767949500d-01,  -2.063496160252500d-02  /)  
            qpta( 4, :) = (/ -2.063496160252500d-02,  -9.587300767949500d-01  /)  
            qpta( 5, :) = (/ -1.258208170141270d-01,  -1.258208170141270d-01  /)  
            qpta( 6, :) = (/ -7.483583659717460d-01,  -1.258208170141270d-01  /)  
            qpta( 7, :) = (/ -1.258208170141270d-01,  -7.483583659717460d-01  /)  
            qpta( 8, :) = (/ -6.235929287619350d-01,  -6.235929287619350d-01  /)  
            qpta( 9, :) = (/ 2.471858575238690d-01,  -6.235929287619350d-01  /)  
            qpta( 10, :) = (/ -6.235929287619350d-01,  2.471858575238690d-01  /)  
            qpta( 11, :) = (/ -9.105409732110950d-01,  -9.105409732110950d-01  /)  
            qpta( 12, :) = (/ 8.210819464221890d-01,  -9.105409732110950d-01  /)  
            qpta( 13, :) = (/ -9.105409732110950d-01,  8.210819464221890d-01  /)  
            qpta( 14, :) = (/ 4.823971975689960d-01,  -5.560740216784690d-01  /)  
            qpta( 15, :) = (/ -5.560740216784690d-01,  4.823971975689960d-01  /)  
            qpta( 16, :) = (/ -9.263231758905270d-01,  -5.560740216784690d-01  /)  
            qpta( 17, :) = (/ -5.560740216784690d-01,  -9.263231758905270d-01  /)  
            qpta( 18, :) = (/ -9.263231758905270d-01,  4.823971975689960d-01  /)  
            qpta( 19, :) = (/ 4.823971975689960d-01,  -9.263231758905270d-01  /)  
            wpta(1) = 1.942715925655980d-01 
            wpta(2) = 6.266940045427800d-02
            wpta(3) = 6.266940045427800d-02
            wpta(4) = 6.266940045427800d-02
            wpta(5) = 1.556550820095490d-01
            wpta(6) = 1.556550820095490d-01
            wpta(7) = 1.556550820095490d-01 
            wpta(8) = 1.592954778544200d-01
            wpta(9) = 1.592954778544200d-01
            wpta(10) = 1.592954778544200d-01
            wpta(11) = 5.115535131739600d-02
            wpta(12) = 5.115535131739600d-02
            wpta(13) = 5.115535131739600d-02
            wpta(14) = 8.656707875457900d-02
            wpta(15) = 8.656707875457900d-02
            wpta(16) = 8.656707875457900d-02
            wpta(17) = 8.656707875457900d-02
            wpta(18) = 8.656707875457900d-02
            wpta(19) = 8.656707875457900d-02

          CASE( 10 )
            nqpta = 25

            qpta( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01  /)  
            qpta( 2, :) = (/ -4.269134091050000d-03,  -4.269134091050000d-03  /)  
            qpta( 3, :) = (/ -9.914617318178990d-01,  -4.269134091050000d-03  /)  
            qpta( 4, :) = (/ -4.269134091050000d-03,  -9.914617318178990d-01  /)  
            qpta( 5, :) = (/ -1.439751005418880d-01,  -1.439751005418880d-01  /)  
            qpta( 6, :) = (/ -7.120497989162250d-01,  -1.439751005418880d-01  /)  
            qpta( 7, :) = (/ -1.439751005418880d-01,  -7.120497989162250d-01  /)  
            qpta( 8, :) = (/ -6.304871745135510d-01,  -6.304871745135510d-01  /)  
            qpta( 9, :) = (/ 2.609743490271020d-01,  -6.304871745135510d-01  /)  
            qpta( 10, :) = (/ -6.304871745135510d-01,  2.609743490271020d-01  /)  
            qpta( 11, :) = (/ -9.590375628566450d-01,  -9.590375628566450d-01  /)  
            qpta( 12, :) = (/ 9.180751257132900d-01,  -9.590375628566450d-01  /)  
            qpta( 13, :) = (/ -9.590375628566450d-01,  9.180751257132900d-01  /)  
            qpta( 14, :) = (/ -7.268528474879330d-01,  6.568468676933890d-01  /)  
            qpta( 15, :) = (/ 6.568468676933890d-01,  -7.268528474879330d-01  /)  
            qpta( 16, :) = (/ -9.299940202054560d-01,  6.568468676933890d-01  /)  
            qpta( 17, :) = (/ 6.568468676933890d-01,  -9.299940202054560d-01  /)  
            qpta( 18, :) = (/ -9.299940202054560d-01,  -7.268528474879330d-01  /)  
            qpta( 19, :) = (/ -7.268528474879330d-01,  -9.299940202054560d-01  /)  
            qpta( 20, :) = (/ -3.345127988227230d-01,  2.594146583058370d-01  /)  
            qpta( 21, :) = (/ 2.594146583058370d-01,  -3.345127988227230d-01  /)  
            qpta( 22, :) = (/ -9.249018594831150d-01,  2.594146583058370d-01  /)  
            qpta( 23, :) = (/ 2.594146583058370d-01,  -9.249018594831150d-01  /)  
            qpta( 24, :) = (/ -9.249018594831150d-01,  -3.345127988227230d-01  /)  
            qpta( 25, :) = (/ -3.345127988227230d-01,  -9.249018594831150d-01  /) 
            wpta(1) = 1.670467996103930d-01
            wpta(2) = 1.445970118411300d-02
            wpta(3) = 1.445970118411300d-02
            wpta(4) = 1.445970118411300d-02
            wpta(5) = 1.489843558419610d-01
            wpta(6) = 1.489843558419610d-01
            wpta(7) = 1.489843558419610d-01 
            wpta(8) = 1.572929468062170d-01
            wpta(9) = 1.572929468062170d-01
            wpta(10) = 1.572929468062170d-01
            wpta(11) = 1.385664617421500d-02
            wpta(12) = 1.385664617421500d-02
            wpta(13) = 1.385664617421500d-02
            wpta(14) = 5.903664066955900d-02
            wpta(15) = 5.903664066955900d-02
            wpta(16) = 5.903664066955900d-02
            wpta(17) = 5.903664066955900d-02
            wpta(18) = 5.903664066955900d-02
            wpta(19) = 5.903664066955900d-02
            wpta(20) = 7.915873439212200d-02
            wpta(21) = 7.915873439212200d-02
            wpta(22) = 7.915873439212200d-02
            wpta(23) = 7.915873439212200d-02
            wpta(24) = 7.915873439212200d-02
            wpta(25) = 7.915873439212200d-02

           CASE DEFAULT

             PRINT*, 'Error: invalid p'
             PRINT*, 'Area quadrature points not available'
             STOP
 
         END SELECT      
         
         ELSE IF (eltype == -1) THEN
         
!            IF (p == 2) THEN
!              qpta(1,:) = (/ 0d0 , (sqrt(3d0)+sqrt(15d0))/6d0 /)
!              qpta(2,:) = (/ 0d0 , (sqrt(3d0)-sqrt(15d0))/6d0 /)
!              qpta(3,:) = (/  sqrt(15d0)/5d0 , (sqrt(87d0)-2d0*sqrt(3d0))/15d0 /)
!              qpta(4,:) = (/ -sqrt(15d0)/5d0 , (sqrt(87d0)-2d0*sqrt(3d0))/15d0 /)
!              qpta(5,:) = (/  sqrt(15d0)/5d0 , (-sqrt(87d0)-2d0*sqrt(3d0))/15d0 /)
!              qpta(6,:) = (/ -sqrt(15d0)/5d0 , (-sqrt(87d0)-2d0*sqrt(3d0))/15d0 /)
!              
!              wpta(1) = 4d0*(2d0/9d0  - (2d0*sqrt(5d0))/45d0)
!              wpta(2) = 4d0*(2d0/9d0  + (2d0*sqrt(5d0))/45d0)
!              wpta(3) = 4d0*(5d0/36d0 + (5d0*sqrt(29d0))/(18d0*29d0))
!              wpta(4) = 4d0*(5d0/36d0 + (5d0*sqrt(29d0))/(18d0*29d0))
!              wpta(5) = 4d0*(5d0/36d0 - (5d0*sqrt(29d0))/(18d0*29d0))
!              wpta(6) = 4d0*(5d0/36d0 - (5d0*sqrt(29d0))/(18d0*29d0))    
!              
!              nqpta = 6
             
!              qpta(1,:) = (/  0.0d0            ,  0.774596669241d0 /)
!              qpta(2,:) = (/  0.563604836881d0 , -0.795508520349d0 /)
!              qpta(3,:) = (/  0.838331011044d0 ,  0.845091361153d0 /)
!              qpta(4,:) = (/  0.651030930900d0 ,  0.166755021097d0 /)
!              qpta(5,:) = (/ -0.484792881050d0 , -0.927694708202d0 /)
!              qpta(6,:) = (/ -0.914603935097d0 , -0.520771886130d0 /)
!              qpta(7,:) = (/ -0.135220856964d0 , -0.279191827433d0 /)
!              qpta(8,:) = (/ -0.731697727745d0 ,  0.417391901524d0 /)           
!              qpta(9,:) = (/ -0.887824220291d0 ,  1.075479856096d0 /)
!              qpta(10,:) = (/ 1.101172842321d0,  -0.485302501018d0 /)
!             
!              
!              wpta(1) = 4d0*0.140845070423d0 
!              wpta(2) = 4d0*0.113931725656d0 
!              wpta(3) = 4d0*0.049023075184d0 
!              wpta(4) = 4d0*0.168918151204d0 
!              wpta(5) = 4d0*0.063463914536d0  
!              wpta(6) = 4d0*0.066611011696d0 
!              wpta(7) = 4d0*0.214897708035d0
!              wpta(8) = 4d0*0.145149421990d0 
!              wpta(9) = 4d0*0.014704280797d0 
!              wpta(10) = 4d0*0.022455640481d0
             
!              nqpta = 10             
!            ELSE
         
             CALL guass_qpts(p,npt,w,r)
           
             n = 1
             DO i = 1,npt
               DO j = 1,npt           
                 wpta(n) = w(i)*w(j) 
                 qpta(n,1) = r(i)
                 qpta(n,2) = r(j)
               
                 n = n + 1
               ENDDO 
             ENDDO
           
             nqpta = npt*npt
           ENDIF
           
!          ENDIF
       
       END SUBROUTINE cubature
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE guass_qpts(p,nqpte,wpte,qpte)

      USE globals, ONLY: pres
      
      IMPLICIT NONE
      INTEGER :: p
      INTEGER :: nqpte
      REAL(pres) :: wpte(13),qpte(13)

      wpte = 0d0
      qpte = 0d0
      
      SELECT CASE(p)

      CASE(2)
         nqpte = 2

         qpte(1) = -0.57735026918963D0
         qpte(2) =  0.57735026918963D0
         
         wpte(1) = 1.D0
         wpte(2) = 1.D0
         
      CASE(3)
         nqpte = 3
         
         qpte(1) = -0.77459666924148D0
         qpte(2) =  0.D0
         qpte(3) =  0.77459666924148D0
         
         wpte(1) =  0.55555555555556D0
         wpte(2) =  0.88888888888888D0
         wpte(3) =  0.55555555555556D0
         
      CASE(4)
         nqpte = 4
         
         qpte(1) = -0.86113631159405D0
         qpte(2) = -0.33998104358486D0
         qpte(3) =  0.33998104358486D0
         qpte(4) =  0.86113631159405D0
         
         wpte(1) =  0.34785484513745D0
         wpte(2) =  0.65214515486255D0
         wpte(3) =  0.65214515486255D0
         wpte(4) =  0.34785484513745D0

      CASE(5)
         nqpte = 5

         qpte(1) = -0.90617984593866D0
         qpte(2) = -0.53846931010568D0
         qpte(3) =  0.D0
         qpte(4) =  0.53846931010568D0
         qpte(5) =  0.90617984593866D0
         
         wpte(1) =  0.23692688505619D0
         wpte(2) =  0.47862867049937D0
         wpte(3) =  0.56888888888889D0
         wpte(4) =  0.47862867049937D0
         wpte(5) =  0.23692688505619D0

      CASE(6)
         nqpte = 6

         qpte(1) = -0.93246951420315D0
         qpte(2) = -0.66120938646626D0
         qpte(3) = -0.23861918608320D0
         qpte(4) =  0.23861918608320D0
         qpte(5) =  0.66120938646626D0
         qpte(6) =  0.93246951420315D0
         
         wpte(1) =  0.17132449237917D0
         wpte(2) =  0.36076157304814D0
         wpte(3) =  0.46791393457269D0
         wpte(4) =  0.46791393457269D0
         wpte(5) =  0.36076157304814D0
         wpte(6) =  0.17132449237917D0

      CASE(7)
         nqpte = 7
         
         qpte(1) = -0.94910791234276D0
         qpte(2) = -0.74153118559939D0
         qpte(3) = -0.40584515137740D0
         qpte(4) =  0.D0
         qpte(5) =  0.40584515137740D0
         qpte(6) =  0.74153118559939D0
         qpte(7) =  0.94910791234276D0
         
         wpte(1) =  0.12948496616887D0
         wpte(2) =  0.27970539148928D0
         wpte(3) =  0.38183005050512D0
         wpte(4) =  0.41795918367347D0
         wpte(5) =  0.38183005050512D0
         wpte(6) =  0.27970539148928D0
         wpte(7) =  0.12948496616887D0

      CASE(8)
         nqpte = 8
         
         qpte(1) = -0.96028985649754D0
         qpte(2) = -0.79666647741363D0
         qpte(3) = -0.52553240991633D0
         qpte(4) = -0.18343464249565D0
         qpte(5) =  0.18343464249565D0
         qpte(6) =  0.52553240991633D0
         qpte(7) =  0.79666647741363D0
         qpte(8) =  0.96028985649754D0
         
         wpte(1) =  0.10122853629038D0
         wpte(2) =  0.22238103445337D0
         wpte(3) =  0.31370664587789D0
         wpte(4) =  0.36268378337836D0
         wpte(5) =  0.36268378337836D0
         wpte(6) =  0.31370664587789D0
         wpte(7) =  0.22238103445337D0
         wpte(8) =  0.10122853629038D0

      CASE(9)
         nqpte = 9
         
         qpte(1) = -0.96816023950763D0
         qpte(2) = -0.83603110732664D0
         qpte(3) = -0.61337143270059D0
         qpte(4) = -0.32425342340381D0
         qpte(5) =  0.D0
         qpte(6) =  0.32425342340381D0
         qpte(7) =  0.61337143270059D0
         qpte(8) =  0.83603110732664D0
         qpte(9) =  0.96816023950763D0
         
         wpte(1) = 0.08127438836163D0
         wpte(2) = 0.18064816069483D0
         wpte(3) = 0.26061069640294D0
         wpte(4) = 0.31234707704000D0
         wpte(5) = 0.33023935500126D0
         wpte(6) = 0.31234707704000D0
         wpte(7) = 0.26061069640294D0
         wpte(8) = 0.18064816069483D0
         wpte(9) = 0.08127438836163D0
         
      CASE(10)
         nqpte = 10
         
         qpte(1)  = -0.97390652851717D0
         qpte(2)  = -0.86506336668898D0
         qpte(3)  = -0.67940956829902D0
         qpte(4)  = -0.43339539412925D0
         qpte(5)  = -0.14887433898163D0
         qpte(6)  =  0.14887433898163D0
         qpte(7)  =  0.43339539412925D0
         qpte(8)  =  0.67940956829902D0
         qpte(9)  =  0.86506336668898D0
         qpte(10) =  0.97390652851717D0

         wpte(1)  =  0.06667134430869D0
         wpte(2)  =  0.14945134915058D0
         wpte(3)  =  0.21908636251598D0
         wpte(4)  =  0.26926671931000D0
         wpte(5)  =  0.29552422471475D0
         wpte(6)  =  0.29552422471475D0
         wpte(7)  =  0.26926671931000D0
         wpte(8)  =  0.21908636251598D0
         wpte(9)  =  0.14945134915058D0
         wpte(10) =  0.06667134430869D0

      CASE(11)
         nqpte = 11
         
         qpte(1)  = -0.97822865814606D0
         qpte(2)  = -0.88706259976810D0
         qpte(3)  = -0.73015200557405D0
         qpte(4)  = -0.51909612920681D0
         qpte(5)  = -0.26954315595234D0
         qpte(6)  =  0.D0
         qpte(7)  =  0.26954315595234D0
         qpte(8)  =  0.51909612920681D0
         qpte(9)  =  0.73015200557405D0
         qpte(10) =  0.88706259976810D0
         qpte(11) =  0.97822865814606D0
         
         wpte(1)  =  0.05566856711627D0
         wpte(2)  =  0.12558036946485D0
         wpte(3)  =  0.18629021092774D0
         wpte(4)  =  0.23319376459199D0
         wpte(5)  =  0.26280454451025D0
         wpte(6)  =  0.27292508677790D0
         wpte(7)  =  0.26280454451025D0
         wpte(8)  =  0.23319376459199D0
         wpte(9)  =  0.18629021092774D0
         wpte(10) =  0.12558036946485D0
         wpte(11) =  0.05566856711627D0
         
      CASE(12)
         nqpte = 12
         
         qpte(1)  = -0.98156063424672D0
         qpte(2)  = -0.90411725637047D0
         qpte(3)  = -0.76990267419430D0
         qpte(4)  = -0.58731795428662D0
         qpte(5)  = -0.36783149899818D0
         qpte(6)  = -0.12523340851147D0
         qpte(7)  =  0.12523340851147D0
         qpte(8)  =  0.36783149899818D0
         qpte(9)  =  0.58731795428662D0
         qpte(10) =  0.76990267419430D0
         qpte(11) =  0.90411725637047D0
         qpte(12) =  0.98156063424672D0
         
         wpte(1)  =  0.04717533638677D0
         wpte(2)  =  0.10693932599520D0
         wpte(3)  =  0.16007832854334D0
         wpte(4)  =  0.20316742672308D0
         wpte(5)  =  0.23349253653835D0
         wpte(6)  =  0.24914704581340D0
         wpte(7)  =  0.24914704581340D0
         wpte(8)  =  0.23349253653835D0
         wpte(9)  =  0.20316742672308D0
         wpte(10) =  0.16007832854334D0
         wpte(11) =  0.10693932599520D0
         wpte(12) =  0.04717533638677D0

      CASE(13)
         nqpte = 13
         
         qpte(1)  = -0.98418305471859D0
         qpte(2)  = -0.91759839922298D0
         qpte(3)  = -0.80157809073331D0
         qpte(4)  = -0.64234933944034D0
         qpte(5)  = -0.44849275103645D0
         qpte(6)  = -0.23045831595513D0
         qpte(7)  =  0.D0
         qpte(8)  =  0.23045831595513D0
         qpte(9)  =  0.44849275103645D0
         qpte(10) =  0.64234933944034D0
         qpte(11) =  0.80157809073331D0
         qpte(12) =  0.91759839922298D0
         qpte(13) =  0.98418305471859D0

         wpte(1)  =  0.04048400476532D0
         wpte(2)  =  0.09212149983773D0
         wpte(3)  =  0.13887351021979D0
         wpte(4)  =  0.17814598076195D0
         wpte(5)  =  0.20781604753689D0
         wpte(6)  =  0.22628318026290D0
         wpte(7)  =  0.23255155323087D0
         wpte(8)  =  0.22628318026290D0
         wpte(9)  =  0.20781604753689D0
         wpte(10) =  0.17814598076195D0
         wpte(11) =  0.13887351021979D0
         wpte(12) =  0.09212149983773D0
         wpte(13) =  0.04048400476532D0

       CASE DEFAULT

         PRINT*, 'Error: invalid p'
         PRINT*, 'Edge quadrature points not available'
         STOP

      END SELECT


      RETURN
      END SUBROUTINE guass_qpts             