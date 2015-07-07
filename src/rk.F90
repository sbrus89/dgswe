      SUBROUTINE rk()

      USE globals, ONLY: dof,ndof,el,ne,nel_type,blk,npart,npartet,elblk, &
                         t,tstage,dt,pt3333,ramp,dramp,ark,brk,crk,rk_type, &
                         Hold,H,MirhsH,Zold,Z,MirhsZ,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy                                            

      IMPLICIT NONE
      
      INTEGER :: et,stg
      
   
!       CALL nan_check()
!       PRINT*, "NaN checked initial condition"      

      rk_type = 22
      SELECT CASE (rk_type)
      
        CASE(11)
        
          CALL swap()
          CALL forward_euler()
          
        CASE(22)
        
          CALL swap()
          CALL forward_euler()

          ! Evaluate RHS
          tstage = t + dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!           CALL rhs3()
!           CALL rhs4()
      
          ! Second RK stage
          DO blk = 1,npart
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN   
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP                
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = .5d0*(Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
                    Z(el,dof) = .5d0*(Zold(el,dof) + Z(el,dof) + dt*MirhsZ(el,dof))
                  ENDDO
!DIR$ IVDEP                  
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = .5d0*(Qxold(el,dof) + Qx(el,dof) + dt*MirhsQx(el,dof))
                  ENDDO
!DIR$ IVDEP                  
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = .5d0*(Qyold(el,dof) + Qy(el,dof) + dt*MirhsQy(el,dof))
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO      

        CASE(33)
        
          CALL swap()
          CALL forward_euler()

          ! Evaluate RHS
          tstage = t + dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!          CALL rhs3()
!          CALL rhs4()
      
          ! Second RK stage
          DO blk = 1,npart
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN    
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = .25d0*(3d0*Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
                    Z(el,dof) = .25d0*(3d0*Zold(el,dof) + Z(el,dof) + dt*MirhsZ(el,dof))
                  ENDDO
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = .25d0*(3d0*Qxold(el,dof) + Qx(el,dof) + dt*MirhsQx(el,dof))
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = .25d0*(3d0*Qyold(el,dof) + Qy(el,dof) + dt*MirhsQy(el,dof))
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO   
      
!          CALL nan_check()      

          ! Evaluate RHS
          tstage = t + .5d0*dt
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
!           CALL rhs3()
!           CALL rhs4()
      
          ! Third RK stage
          DO blk = 1,npart
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN      
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP      
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
!                     H(el,dof) = pt3333*(Hold(el,dof) + 2d0*(H(el,dof) + dt*MirhsH(el,dof)))
                    Z(el,dof) = pt3333*(Zold(el,dof) + 2d0*(Z(el,dof) + dt*MirhsZ(el,dof)))
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qx(el,dof) = pt3333*(Qxold(el,dof) + 2d0*(Qx(el,dof) + dt*MirhsQx(el,dof)))
                  ENDDO
!DIR$ IVDEP        
!!DIR$ VECTOR ALIGNED
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qy(el,dof) = pt3333*(Qyold(el,dof) + 2d0*(Qy(el,dof) + dt*MirhsQy(el,dof)))
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO      

        CASE(45)
                    
          DO stg = 1,5

          tstage = t + crk(stg)*dt 
          ramp = TANH((2d0*tstage)/(86400d0*dramp))
          CALL rhs2()
      
          DO blk = 1,npart
            DO et  = 1,nel_type
              IF (npartet(et,blk) > 0) THEN    
          
                DO dof = 1,ndof(et)
!DIR$ IVDEP
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Zold(el,dof) = ark(stg)*Zold(el,dof) + dt*MirhsZ(el,dof)
                    Z(el,dof) = Z(el,dof) + brk(stg)*Zold(el,dof)
                  ENDDO
!DIR$ IVDEP
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qxold(el,dof) = ark(stg)*Qxold(el,dof) + dt*MirhsQx(el,dof)
                    Qx(el,dof) = Qx(el,dof) + brk(stg)*Qxold(el,dof)
                  ENDDO
!DIR$ IVDEP        
                  DO el = elblk(1,blk,et),elblk(2,blk,et)
                    Qyold(el,dof) = ark(stg)*Qyold(el,dof) + dt*MirhsQy(el,dof)
                    Qy(el,dof) = Qy(el,dof) + brk(stg)*Qyold(el,dof)
                  ENDDO
                ENDDO
      
              ENDIF
            ENDDO
          ENDDO   
      
          ENDDO
          
        CASE DEFAULT
        
          PRINT*, "Time-stepping option not availiable: ", rk_type
          
      END SELECT          
      
      

      CALL nan_check()
         

      RETURN
      END SUBROUTINE RK
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE swap()
      
      USE globals, ONLY: dof,ndof,el,nel_type,blk,npart,npartet,elblk, &
                         Hold,H,Zold,Z,Qxold,Qx,Qyold,Qy      
                       
      IMPLICIT NONE                       
                       
      INTEGER :: et                         
                         
!       CALL nan_check()
!       PRINT*, "NaN checked initial condition"      

      ! Save previous solution
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP            
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 Hold(el,dof) = H(el,dof)
                Zold(el,dof) = Z(el,dof)
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxold(el,dof) = Qx(el,dof)
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qyold(el,dof) = Qy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO                         
      
      END SUBROUTINE swap
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE forward_euler()
      
      USE globals, ONLY: dof,ndof,el,ne,nel_type,blk,npart,npartet,elblk, &
                         t,tstage,dt,ramp,dramp, &
                         Hold,H,MirhsH,Zold,Z,MirhsZ,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy        
      
      IMPLICIT NONE
      
      INTEGER :: et
      
      ! Evalutate right hand side  
      tstage = t
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
!       CALL rhs3()
!        CALL rhs4()
      
      ! First RK stage
      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP            
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 H(el,dof) = Hold(el,dof) + dt*MirhsH(el,dof)
                Z(el,dof) = Zold(el,dof) + dt*MirhsZ(el,dof)
              ENDDO
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = Qxold(el,dof) + dt*MirhsQx(el,dof)
              ENDDO
!DIR$ IVDEP              
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = Qyold(el,dof) + dt*MirhsQy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO
           
      
      END SUBROUTINE forward_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE nan_check()
      
      USE globals, ONLY: dof,ndof,el,ne, &
                         H,Z,Qx,Qy, &
                         blk,npart,nel_type,elblk,npartet, &
                         t
                         
      USE messenger2, ONLY: finish                               
      
      IMPLICIT NONE
      
      INTEGER :: et
      
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
      
            DO dof = 1,ndof(et)
      
              DO el = elblk(1,blk,et),elblk(2,blk,et)
!                 IF (H(el,dof) /= H(el,dof)) THEN
                IF (Z(el,dof) /= Z(el,dof)) THEN
                  PRINT*, "NaN detected in H solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  CALL write_output()
                  CALL abort()
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qx(el,dof) /= Qx(el,dof)) THEN
                  PRINT*, "NaN detected in Qx solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  CALL write_output()
                  CALL abort()
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qy(el,dof) /= Qy(el,dof)) THEN
                  PRINT*, "NaN detected in Qy solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  CALL write_output()
                  CALL abort()
                ENDIF
              ENDDO
        
            ENDDO
      
          ENDIF
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE nan_check
