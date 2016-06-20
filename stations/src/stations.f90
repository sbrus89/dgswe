      PROGRAM stations

      USE globals, ONLY:rp,ndof,mndof,np,mnnds,nnds,nverts,nel_type,ne,nn,nepn,mnepn,epn, &                                    
                        xy,ect,el_type,elxy,elhb,H,Qx,Qy,elsta,hbsta,nsta,xysta,phi_sta,t
      USE read_dginp, ONLY:read_input,p,ctp,hbp,out_direc,sol_opt,sol_snap,sta_opt
      USE edge_connectivity_mod, ONLY: elements_per_node
      USE messenger2, ONLY: message_init


      

      IMPLICIT NONE

      INTEGER :: el,sta,line,lines,dof
      INTEGER :: elin,et
      REAL(rp) :: Hsta,Qxsta,Qysta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: l,phi    
      CHARACTER(100) :: tmp
      


      CALL message_init()
      CALL read_input()
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)*(p+1)
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)
      mndof = maxval(ndof)
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp        
      
      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      mnnds = maxval(nnds)     
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4          
                        
      lines = int(sol_snap)                  
                        
      CALL read_grid()                 
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)                   

      
      CALL create_101btype_stations()
      
      CALL find_stations()
      
      STOP      
      
        
        
        
        

      
      ALLOCATE(H(ne,mndof),Qx(ne,mndof),Qy(ne,mndof))   

      
      OPEN(UNIT=611,FILE=trim(out_direc) //"station_H.d")
      OPEN(UNIT=612,FILE=trim(out_direc) //"station_hb.d")      
      OPEN(UNIT=621,FILE=trim(out_direc) //"station_Qx.d")     
      OPEN(UNIT=622,FILE=trim(out_direc) //"station_Qy.d")           
      
      OPEN(UNIT=63, FILE=trim(out_direc) //"solution_H.d")
      OPEN(UNIT=641,FILE=trim(out_direc) //"solution_Qx.d")     
      OPEN(UNIT=642,FILE=trim(out_direc) //"solution_Qy.d")      
      
      READ(63,*) tmp
      READ(641,*) tmp 
      READ(642,*) tmp
      
      WRITE(611,"(A)") tmp
      WRITE(612,"(A)") tmp      
      WRITE(621,"(A)") tmp
      WRITE(622,"(A)") tmp      
      
      WRITE(611,*) nsta
      WRITE(612,*) nsta      
      WRITE(621,*) nsta
      WRITE(622,*) nsta
      
      DO sta = 1,nsta
        WRITE(612,"(E24.17)") hbsta(sta)
      ENDDO      
      
      DO line = 1,lines+1
        READ(63,*) t
        READ(641,*) t
        READ(642,*) t
        
        PRINT("(A,I5)"), "Writing time snap: ", line
        
        DO dof = 1,mndof
          READ(63,*) (H(el,dof), el = 1,ne)
          READ(641,*) (Qx(el,dof), el = 1,ne)
          READ(642,*) (Qy(el,dof), el = 1,ne)
        ENDDO
        
        WRITE(611,*) t
        WRITE(621,*) t
        WRITE(622,*) t           
                
        DO sta = 1,nsta
                 
          elin = elsta(sta)                 
          et = el_type(elin)                 
                    
          Hsta = 0d0
          Qxsta = 0d0
          Qysta = 0d0
          DO dof = 1,ndof(et)
            Hsta = Hsta + H(elin,dof)*phi_sta(dof,sta)
            Qxsta = Qxsta + Qx(elin,dof)*phi_sta(dof,sta)            
            Qysta = Qysta + Qy(elin,dof)*phi_sta(dof,sta)
          ENDDO                
          
          WRITE(611,"(E24.17)") Hsta
          WRITE(621,"(E24.17)") Qxsta
          WRITE(622,"(E24.17)") Qysta
        ENDDO
        
      ENDDO
      
      CLOSE(63)
      CLOSE(641)
      CLOSE(642)
      
      CLOSE(611)
      CLOSE(621)
      CLOSE(622)      
      END PROGRAM stations