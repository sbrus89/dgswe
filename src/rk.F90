      SUBROUTINE rk()

      USE globals, ONLY: dof,ndof,el,ne, &
                         t,tstage,dt,pt3333,ramp,dramp, &
                         Hold,H,rhsH,Qxold,Qx,rhsQx,Qyold,Qy,rhsQy

      IMPLICIT NONE

      ! Save previous solution
      DO dof = 1,ndof
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Hold(el,dof) = H(el,dof)
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qxold(el,dof) = Qx(el,dof)
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne 
          Qyold(el,dof) = Qy(el,dof)
        ENDDO
      ENDDO

      ! Evalutate right hand side  
      tstage = t
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
      
      ! First RK stage
      DO dof = 1,ndof
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          H(el,dof) = Hold(el,dof) + dt*rhsH(el,dof)
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qx(el,dof) = Qxold(el,dof) + dt*rhsQx(el,dof)
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qy(el,dof) = Qyold(el,dof) + dt*rhsQy(el,dof)
        ENDDO
      ENDDO

#ifdef rk22
      ! Evaluate RHS
      tstage = t + dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
      
      ! Second RK stage
      DO dof = 1,ndof
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          H(el,dof) = .5d0*(Hold(el,dof) + H(el,dof) + dt*rhsH(el,dof))
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qx(el,dof) = .5d0*(Qxold(el,dof) + Qx(el,dof) + dt*rhsQx(el,dof))
        ENDDO
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qy(el,dof) = .5d0*(Qyold(el,dof) + Qy(el,dof) + dt*rhsQy(el,dof))
        ENDDO
      ENDDO
#endif



#ifdef rk33
! 3rd order TVD-RK (1st stage is the same as 2nd order)

      ! Evaluate RHS
      tstage = t + dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
      
      ! Second RK stage
      DO dof = 1,ndof
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          H(el,dof) = .25d0*(3d0*Hold(el,dof) + H(el,dof) + dt*rhsH(el,dof))
        ENDDO
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qx(el,dof) = .25d0*(3d0*Qxold(el,dof) + Qx(el,dof) + dt*rhsQx(el,dof))
        ENDDO
!DIR$ IVDEP        
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qy(el,dof) = .25d0*(3d0*Qyold(el,dof) + Qy(el,dof) + dt*rhsQy(el,dof))
        ENDDO
      ENDDO

      ! Evaluate RHS
      tstage = t + .5d0*dt
      ramp = TANH((2d0*tstage)/(86400d0*dramp))
      CALL rhs2()
      
      ! Third RK stage
      DO dof = 1,ndof
!DIR$ IVDEP      
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          H(el,dof) = pt3333*(Hold(el,dof) + 2d0*(H(el,dof) + dt*rhsH(el,dof)))
        ENDDO
!DIR$ IVDEP        
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qx(el,dof) = pt3333*(Qxold(el,dof) + 2d0*(Qx(el,dof) + dt*rhsQx(el,dof)))
        ENDDO
!DIR$ IVDEP        
!DIR$ VECTOR ALIGNED
        DO el = 1,ne
          Qy(el,dof) = pt3333*(Qyold(el,dof) + 2d0*(Qy(el,dof) + dt*rhsQy(el,dof)))
        ENDDO
      ENDDO
#endif

!       CALL nan_check()


      RETURN
      END SUBROUTINE RK
      
      
      
      
      SUBROUTINE nan_check()
      
      USE globals, ONLY: dof,ndof,el,ne, &
                         H,Qx,Qy
      
      IMPLICIT NONE
      
      DO dof = 1,ndof
      
        DO el = 1,ne
          IF (H(el,dof) /= H(el,dof)) THEN
            PRINT*, "NaN detected in H solution"
            STOP
          ENDIF
        ENDDO
        
        DO el = 1,ne
          IF (Qx(el,dof) /= Qx(el,dof)) THEN
            PRINT*, "NaN detected in Qx solution"
            STOP
          ENDIF
        ENDDO
        
        DO el = 1,ne
          IF (Qy(el,dof) /= Qy(el,dof)) THEN
            PRINT*, "NaN detected in Qy solution"
            STOP
          ENDIF
        ENDDO
        
      ENDDO
      
      RETURN
      END SUBROUTINE nan_check