      SUBROUTINE rk()

      USE globals, ONLY: dof,ndof,el,ne,nel_type,blk,npart,npartet,elblk, &
                         t,tstage,dt,pt3333,ramp,dramp, &
                         Hold,H,MirhsH,Qxold,Qx,MirhsQx,Qyold,Qy,MirhsQy

      IMPLICIT NONE
      
      INTEGER :: et
      

!       CALL nan_check()
!       PRINT*, "NaN checked initial condition"      

      ! Save previous solution
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
          
            DO dof = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Hold(el,dof) = H(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qxold(el,dof) = Qx(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qyold(el,dof) = Qy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO

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
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = Hold(el,dof) + dt*MirhsH(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = Qxold(el,dof) + dt*MirhsQx(el,dof)
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = Qyold(el,dof) + dt*MirhsQy(el,dof)
              ENDDO
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO

#ifdef rk22
      ! Evaluate RHS
      tstage = t + dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
!       CALL rhs3()
!        CALL rhs4()
      
      ! Second RK stage
      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN   
          
            DO dof = 1,ndof(et)
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = .5d0*(Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qx(el,dof) = .5d0*(Qxold(el,dof) + Qx(el,dof) + dt*MirhsQx(el,dof))
              ENDDO
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                Qy(el,dof) = .5d0*(Qyold(el,dof) + Qy(el,dof) + dt*MirhsQy(el,dof))
              ENDDO
            ENDDO
      
          ENDIF
        ENDDO
      ENDDO      
#endif



#ifdef rk33
! 3rd order TVD-RK (1st stage is the same as 2nd order)

      ! Evaluate RHS
      tstage = t + dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
       CALL rhs2()
!      CALL rhs3()
!        CALL rhs4()
      
      ! Second RK stage
      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN    
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = .25d0*(3d0*Hold(el,dof) + H(el,dof) + dt*MirhsH(el,dof))
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
      
      

      ! Evaluate RHS
      tstage = t + .5d0*dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
!       CALL rhs3()
!       CALL rhs4()
      
      ! Third RK stage
      DO blk = 1,npart
        DO et  = 1,nel_type
          IF (npartet(et,blk) > 0) THEN      
          
            DO dof = 1,ndof(et)
!DIR$ IVDEP      
!!DIR$ VECTOR ALIGNED
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                H(el,dof) = pt3333*(Hold(el,dof) + 2d0*(H(el,dof) + dt*MirhsH(el,dof)))
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
#endif

     CALL nan_check()

      RETURN
      END SUBROUTINE RK
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE nan_check()
      
      USE globals, ONLY: dof,ndof,el,ne, &
                         H,Qx,Qy, &
                         blk,npart,nel_type,elblk,npartet, &
                         t
      
      IMPLICIT NONE
      
      INTEGER :: et
      
      DO blk = 1,npart
        DO et = 1,nel_type
          IF (npartet(et,blk) > 0) THEN
      
            DO dof = 1,ndof(et)
      
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (H(el,dof) /= H(el,dof)) THEN
                  PRINT*, "NaN detected in H solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  STOP
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qx(el,dof) /= Qx(el,dof)) THEN
                  PRINT*, "NaN detected in Qx solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  STOP
                ENDIF
              ENDDO
        
              DO el = elblk(1,blk,et),elblk(2,blk,et)
                IF (Qy(el,dof) /= Qy(el,dof)) THEN
                  PRINT*, "NaN detected in Qy solution"
                  PRINT("(A,e15.8)"), 't = ', t
                  STOP
                ENDIF
              ENDDO
        
            ENDDO
      
          ENDIF
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE nan_check
