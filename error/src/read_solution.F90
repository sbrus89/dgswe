      SUBROUTINE read_solution(sol,t)

      USE globals, ONLY:rp,solution,lines
      USE read_write_output, ONLY: read_solution_full,read_fort6163,read_fort6264

      IMPLICIT NONE
      
      TYPE(solution) :: sol
      
      INTEGER :: i,dof,el,line,et,n,nd
      INTEGER :: nsnap_read
      REAL(rp) :: t,hb
      REAL(rp), DIMENSION(:), ALLOCATABLE :: t_sol
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: eta,u,v
      
      PRINT("(A)"), TRIM(sol%sol_name) // " solution"
      
#ifdef dgswem      
      
      OPEN(UNIT=11, FILE=trim(sol%out_direc) //"DG.63")
      OPEN(UNIT=12, FILE=trim(sol%out_direc) //"DG.64")     
      
      READ(11,*) 
      READ(11,*)
      READ(12,*)
      READ(12,*)
      
      DO line = 1,lines
        READ(11,*) t,n
        READ(12,*) t,n 
        
        DO el = 1,sol%ne
          et = sol%el_type(el)
          DO dof = 1,sol%ndof(et)
            READ(11,*) i,sol%H(el,dof) !, hb
            READ(12,*) i,sol%Qx(el,dof), sol%Qy(el,dof)
          ENDDO    
        ENDDO        
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12) 
      
#elif adcirc     

      nsnap_read = lines
      CALL read_fort6163(sol%out_direc,"63",t_sol,eta,nsnap_read,last_snap="T")
      CALL read_fort6264(sol%out_direc,"64",t_sol,u,v,nsnap_read,last_snap="T")
      
      ALLOCATE(sol%H(sol%ne,3,1),sol%Qx(sol%ne,3,1),sol%Qy(sol%ne,3,1))
      DO el = 1,sol%ne
        DO nd = 1,3
          sol%H(el,nd,1) = eta(sol%ect(nd,el),1)
          sol%Qx(el,nd,1) = u(sol%ect(nd,el),1)
          sol%Qy(el,nd,1) = v(sol%ect(nd,el),1)
        ENDDO
      ENDDO
      
#else      
      
      
!       OPEN(UNIT=11, FILE=trim(sol%out_direc) //"solution_H.d")
!       OPEN(UNIT=12, FILE=trim(sol%out_direc) //"solution_Qx.d")     
!       OPEN(UNIT=13, FILE=trim(sol%out_direc) //"solution_Qy.d")      
!       
!       READ(11,*) 
!       READ(12,*)
!       READ(13,*)
!       
!       DO line = 1,lines+1
!         READ(11,*) t
!         READ(12,*) t
!         READ(13,*) t 
!         
!         DO dof = 1,sol%mndof
!           READ(11,*) (sol%H(el,dof), el = 1,sol%ne)
!           READ(12,*) (sol%Qx(el,dof), el = 1,sol%ne)
!           READ(13,*) (sol%Qy(el,dof), el = 1,sol%ne)
!         ENDDO    
!         
!       ENDDO
!       
!       CLOSE(11)
!       CLOSE(12)
!       CLOSE(13)     
      
      nsnap_read = lines+1
      CALL read_solution_full(sol%out_direc,'Z.sol',"T",t_sol,sol%H,nsnap_read,last_snap="T")       
      CALL read_solution_full(sol%out_direc,'Qx.sol',"T",t_sol,sol%Qx,nsnap_read,last_snap="T")           
      CALL read_solution_full(sol%out_direc,'Qy.sol',"T",t_sol,sol%Qy,nsnap_read,last_snap="T")           

#endif      

      PRINT("(A)"), "Done"
      PRINT*, " "

      RETURN 
      END SUBROUTINE read_solution
